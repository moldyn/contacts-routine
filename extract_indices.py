import click
import numpy as np

from txt_io import savetxt


@click.command(
    no_args_is_help='-h',
    help='Plot average formed contacts',
)
@click.option(
    '--is-contacts',
    'iscontactfile',
    required=True,
    type=click.Path(exists=True),
    help='Path to contacts file',
)
@click.option(
    '--threshold',
    required=True,
    type=click.FloatRange(min=0, max=1),
    help='Threshold of fraction of formed contacts to be selected.',
)
@click.option(
    '--max-threshold',
    type=click.FloatRange(min=0, max=1),
    default=1.0,
    help='Maximum threshold of fraction of formed contacts (default: 1.0 = no upper limit).',
)
def main(iscontactfile, threshold, max_threshold):
    # Validate thresholds
    if threshold > max_threshold:
        raise click.UsageError(
            f"Minimum threshold ({threshold}) cannot be greater than maximum threshold ({max_threshold})"
        )
    # load files
    is_contacts = np.loadtxt(iscontactfile, usecols=2, ndmin=1)
    idxs = np.loadtxt(iscontactfile, usecols=(0, 1), dtype=int, ndmin=2)  # starting from 1

    # draw rectangles to highlight native contacts
    selected_idxs = [
        (i, j)
        for (i, j), is_contact in zip(idxs, is_contacts)
        if threshold <= is_contact <= max_threshold
    ]

    output_filename = f'{iscontactfile}.thr{threshold}-{max_threshold}.ndx'
    #if max_threshold < 1.0:
        #output_filename += f'-{max_threshold:g}'
    #output_filename += '.ndx'
    
    savetxt(
        output_filename,
        selected_idxs,
        fmt='%.0f',
    )


if __name__ == '__main__':
    main()
