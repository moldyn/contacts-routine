import itertools

import click
import numpy as np
import mdtraj as md
from tqdm import tqdm

from trajectory_files import TRAJ_EXTENSIONS, get_trajectory_files


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
    '--topology',
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
def main(trajfile, trajlist, topfile, ndxfile, output, trajext):
    # Get list of trajectory files
    traj_files = get_trajectory_files(trajfile, trajlist, topfile, trajext)
    
    print(f"Processing {len(traj_files)} trajectory file(s)")
    if len(traj_files) <= 10:
        for f in traj_files:
            print(f"  - {f}")
    else:
        print(f"  - {traj_files[0]}")
        print(f"  - ...")
        print(f"  - {traj_files[-1]}")
    
    # load files
    index_pairs = np.loadtxt(ndxfile, dtype=int)
    # Validate residue numbering
    top = md.load(topfile).topology
    pdb_resids = sorted(set([res.resSeq for res in top.residues]))

    # Check for negative residues
    if any(r < 0 for r in pdb_resids):
        raise click.UsageError(
            f"PDB contains negative residue indices. Please renumber starting from 1."
        )

    # Check if indices are sequential (required for mdtraj)
    if pdb_resids != list(range(min(pdb_resids), max(pdb_resids) + 1)):
        raise click.UsageError(
            "PDB has non-sequential residue numbering. This script requires "
            "sequential numbering (e.g., 1,2,3,... with no gaps). "
            "Please renumber your PDB, e.g. with: gmx editconf -f in.pdb -o out.pdb -resnr 1"
        )


    # convert residue indices to heavy atoms
    resseq_to_mdtraj = {res.resSeq: res.index for res in top.residues}

    atoms_per_res = {
        index: [
            atom.index for atom in top.residue(resseq_to_mdtraj[index]).atoms
            if atom.element.symbol != 'H'
        ]
        for index in np.unique(index_pairs)
    }

    atom_pairs = [
        list(
            itertools.product(
                atoms_per_res[i], atoms_per_res[j],
            ),
        ) for i, j in index_pairs
    ]

    index_output = np.concatenate([
    [
        (res_i, res_j, atom_i + 1, atom_j + 1)  # only atoms get +1
        for atom_i, atom_j in atom_pair
    ]
    for (res_i, res_j), atom_pair in zip(index_pairs, atom_pairs)
    ])
    np.savetxt(
        f'{output}.atom_indices',
        index_output,  # no global +1
        fmt='%.0f',
        header='res_i res_j atom_i atom_j',
    )

    # flatten atom_pairs to loop over all
    atom_pairs = np.concatenate(atom_pairs)

    with open(output, 'w') as ostream:
        for traj_idx, traj_file in enumerate(traj_files):
            print(f"Processing trajectory {traj_idx + 1}/{len(traj_files)}: {traj_file}")
            for distances in compute_distances(
                traj_file, topfile, atom_pairs,
            ):
                np.savetxt(ostream, distances, fmt='%.5f')


def compute_distances(trajfile, topfile, atom_pairs, chunk=100):
    # iterload reads any mdtraj-supported format chunk by chunk
    for frames in tqdm(
        md.iterload(trajfile, top=topfile, chunk=chunk),
        desc=f"Processing {trajfile}",
    ):
        yield md.compute_distances(
            frames,
            atom_pairs=atom_pairs,
        )


if __name__ == '__main__':
    main()
