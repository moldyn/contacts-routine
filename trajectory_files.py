"""Shared helpers to collect trajectory files of any supported format."""
import glob
import os

import click

# Trajectory formats readable by both MDAnalysis (estimate_contacts.py) and
# mdtraj (contacts.py). The format is detected from the file extension.
TRAJ_EXTENSIONS = (
    'xtc', 'trr', 'dcd', 'nc', 'ncdf', 'mdcrd', 'xyz', 'gro', 'pdb',
)


def find_trajectories(folder, topfile, extension=None):
    """Return sorted trajectory files in folder.

    If extension is None, all files with a supported trajectory extension are
    collected (the topology file itself is skipped). If more than one format
    is found an error is raised, to avoid counting the same simulation twice
    (e.g. traj.xtc and traj.trr).
    """
    extensions = (
        [extension.lstrip('.').lower()] if extension else TRAJ_EXTENSIONS
    )
    topfile = os.path.realpath(topfile)
    traj_files = sorted(
        path
        for path in glob.glob(os.path.join(folder, '*'))
        if os.path.isfile(path)
        and path.rsplit('.', 1)[-1].lower() in extensions
        and os.path.realpath(path) != topfile
    )
    if not traj_files:
        raise click.UsageError(
            f'No trajectory files ({", ".join(extensions)}) found in folder: '
            f'{folder}'
        )

    found = sorted({path.rsplit('.', 1)[-1].lower() for path in traj_files})
    if len(found) > 1:
        raise click.UsageError(
            f'Found trajectory files of multiple formats ({", ".join(found)}) '
            f'in folder: {folder}. Please select one with --traj-ext.'
        )
    return traj_files


def get_trajectory_files(trajfile, trajlist, topfile, extension=None):
    """Resolve --trajectory / --trajectory-list into a list of files."""
    if trajlist:
        if os.path.isdir(trajlist):
            return find_trajectories(trajlist, topfile, extension)
        # It's a file with list of trajectories
        with open(trajlist, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    if trajfile:
        return [trajfile]
    raise click.UsageError(
        'Either --trajectory or --trajectory-list must be provided'
    )
