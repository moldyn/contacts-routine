from functools import partial
from multiprocessing import cpu_count, Pool, Lock

import click
import MDAnalysis
import numpy as np
from MDAnalysis.lib.distances import capped_distance, distance_array
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

from trajectory_files import TRAJ_EXTENSIONS, get_trajectory_files

CUTOFF = 4.5  # [AA]


@click.command(
    no_args_is_help='-h',
    help='Mindist computes the minimal distance between pairs of residues.',
)
@click.option(
    '--trajectory',
    '-f',
    'trajfile',
    type=click.Path(exists=True),
    help='Path to trajectory file (e.g. .xtc, .trr, .dcd, .nc)',
)
@click.option(
    '--trajectory-list',
    '--traj-list',
    'trajlist',
    type=click.Path(exists=True),
    help=(
        'Path to folder containing trajectory files '
        f'({", ".join(TRAJ_EXTENSIONS)}) OR file with list of trajectory paths'
    ),
)
@click.option(
    '--traj-ext',
    'trajext',
    help=(
        'Only use trajectory files with this extension when '
        '--trajectory-list is a folder (e.g. xtc, trr, dcd).'
    ),
)
@click.option(
    '--top',
    '-s',
    'topfile',
    required=True,
    type=click.Path(exists=True),
    help='Path to topology file (.tpr or .pdb)',
)
@click.option(
    '--index',
    '-n',
    'ndxfile',
    required=True,
    type=click.Path(exists=True),
    help='Path to index file, 1 (!) indexed of shape (n, 2).',
)
@click.option(
    '--output',
    '-o',
    required=True,
    type=click.Path(),
    help='Path to output file',
)
@click.option(
    '--count-hydrogen',
    is_flag=True,
    help='Count hydrogen atoms.',
)
@click.option(
    '--mode',
    type=click.Choice(['overall', 'per-trajectory'], case_sensitive=False),
    default='overall',
    help=(
        'Contact frequency calculation mode. '
        'overall: average across all trajectories (default). '
        'per-trajectory: contacts must meet threshold in each trajectory individually.'
    ),
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Print detailed progress information',
)
def main(
    trajfile, trajlist, trajext, topfile, ndxfile, output, count_hydrogen,
    mode, verbose,
):
    # Get list of trajectory files
    traj_files = get_trajectory_files(trajfile, trajlist, topfile, trajext)
    if verbose:
        print(f"Processing {len(traj_files)} trajectory file(s)")
        if len(traj_files) <= 10:
            for f in traj_files:
                print(f"  - {f}")
        else:
            print(f"  - {traj_files[0]}")
            print(f"  - ...")
            print(f"  - {traj_files[-1]}")
    
    contact_pairs = np.loadtxt(ndxfile, dtype=int)
    n_contact_pairs = len(contact_pairs)
    # Validate residue numbering in PDB
    universe_check = MDAnalysis.Universe(topfile, traj_files[0])
    pdb_resids = sorted(set([res.resid for res in universe_check.residues]))

    # Check for negative residues
    if any(r < 0 for r in pdb_resids):
        raise click.UsageError(
            f"PDB contains negative residue indices: {[r for r in pdb_resids if r < 0]}. "
            "Please renumber residues starting from 1."
            )

    # Check for gaps
    expected_resids = set(range(min(pdb_resids), max(pdb_resids) + 1))
    missing = expected_resids - set(pdb_resids)
    if missing:
        click.echo(f"WARNING: PDB has gaps in residue numbering. Missing resids: {sorted(missing)}", err=True)

    # Check if requested indices exist in PDB
    requested_resids = set(contact_pairs.flatten())
    invalid = requested_resids - set(pdb_resids)
    if invalid:
        raise click.UsageError(
            f"Index file contains residue IDs not found in PDB: {sorted(invalid)}. "
            f"PDB residue range: {min(pdb_resids)}-{max(pdb_resids)}"
        )
    select_atoms_str = (
        'resid {res}' if count_hydrogen else 'resid {res} and not type H'
    )

    # Initialize atoms_by_res using first trajectory
    universe = MDAnalysis.Universe(topfile, traj_files[0])
    atoms_by_res = {}
    residues = np.unique(contact_pairs)
    for res in residues:
        atoms_by_res[res] = universe.select_atoms(
            select_atoms_str.format(res=res),
        )

    # Process each trajectory
    total_frames = 0
    cmap_total = np.zeros(len(contact_pairs))
    
    if mode == 'per-trajectory':
        # For per-trajectory mode, track minimum fraction across trajectories
        cmap_per_traj = []
    
    for traj_idx, traj_file in enumerate(traj_files):
        if verbose:    
            print(f"\nProcessing trajectory {traj_idx + 1}/{len(traj_files)}: {traj_file}")
        universe = MDAnalysis.Universe(topfile, traj_file)
        # Update atom selections for this trajectory
        for res in residues:
            atoms_by_res[res] = universe.select_atoms(
                select_atoms_str.format(res=res),
                )
        
        # Define slices for parallel processing
        n_jobs = cpu_count()
        n_frames = universe.trajectory.n_frames
        n_blocks = n_jobs
        n_frames_per_slice = n_frames // n_blocks

        slices_values = [
            range(
                i * n_frames_per_slice,
                (i + 1) * n_frames_per_slice,
            )
            for i in range(n_blocks - 1)
        ]
        slices_values.append(
            range(
                (n_blocks - 1) * n_frames_per_slice,
                n_frames,
            )
        )

        # Loop through trajectory
        run_per_slice = partial(
            cmap_per_traj_slice,
            blockslices=slices_values,
            topfile=topfile,
            trajfile=traj_file,
            contact_pairs=contact_pairs,
            select_atoms_str=select_atoms_str,
        )

        with Pool(
            n_jobs, initializer=tqdm.set_lock, initargs=(Lock(), ),
        ) as workers:
            result = workers.map(run_per_slice, np.arange(n_blocks))

        cmap_traj = np.sum(result, axis=0)
        
        if mode == 'per-trajectory':
            # Store fraction for this trajectory
            cmap_per_traj.append(cmap_traj / n_frames)
        
        cmap_total += cmap_traj
        total_frames += n_frames
        if verbose:
            print(f"Trajectory {traj_idx + 1}: {n_frames} frames processed")

    # Calculate final contact map based on mode
    if mode == 'overall':
        cmap = cmap_total / total_frames
        mode_description = f'averaged over {total_frames} frames from {len(traj_files)} trajectory files'
    else:  # per-trajectory
        # Take minimum fraction across all trajectories
        cmap = np.min(cmap_per_traj, axis=0)
        mode_description = f'minimum across {len(traj_files)} trajectories (must pass threshold in each)'
    if verbose:
        print(f"\nTotal frames processed: {total_frames}")
        print(f"Mode: {mode}")

    np.savetxt(
        output,
        [
            f'{i:.0f} {j:.0f} {cmap[idx]:.5f}'
            for idx, (i, j) in enumerate(contact_pairs)
        ],
        header=(
            f'FC = idx_i idx_j FC ({mode_description})'
        ),
        fmt='%s',
    )


def cmap_per_traj_slice(
    slice_idx, blockslices, topfile, trajfile, contact_pairs, select_atoms_str,
):
    # Each worker creates its own Universe - avoids multiprocessing issues
    universe = MDAnalysis.Universe(topfile, trajfile)
    
    # Create atom selections for this worker
    atoms_by_res = {}
    residues = np.unique(contact_pairs)
    for res in residues:
        atoms_by_res[res] = universe.select_atoms(
            select_atoms_str.format(res=res),
        )
    
    cmap = np.zeros(len(contact_pairs))
    blockslice = blockslices[slice_idx]
    for ts in tqdm(
        universe.trajectory[blockslice.start:blockslice.stop],
        position=slice_idx,
        desc=f'process {slice_idx:>2.0f}',
        leave=False,
        mininterval=1,
    ):
        # read the box of the current frame (it changes in NPT simulations;
        # some readers, e.g. DCD, do not update it in place)
        box = ts.dimensions
        cmap += np.array([
            np.min(
                distance_array(
                    atoms_by_res[resi].positions,
                    atoms_by_res[resj].positions,
                    box=box,
                    backend='serial',
                ),
            ) <= CUTOFF
            for idx, (resi, resj) in enumerate(contact_pairs)
        ]).astype(int)
    return cmap


if __name__ == "__main__":
    main()
