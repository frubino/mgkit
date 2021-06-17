"""
.. versionadded: 0.5.7

"""

import logging
from os import sep
import click
from tqdm import tqdm
from .. import logger
from . import utils
from mgkit.io import gff
from mgkit import align


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('depth', help="""
Reads a [Depth file] made with `samtools depth` and a [GFF file] to output file with information about
annotations mean coverage
""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-g', '--gff-file', default=None, type=click.File('rb'), required=True,
                help="GFF with annotation to search")
@click.option('-s', '--separator', default='\t', type=click.STRING, show_default=True,
              help="Field separator for output file")
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('depth-file', type=click.File('rb', lazy=False), default='-')
@click.argument('output-file', type=click.File('w', lazy=False), default='-')
def depth_command(verbose, gff_file, separator, progress, depth_file, output_file):
    """
    .. versionadded:: 0.5.7
    """

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    max_cov = 0
    annotations = {}

    for annotation in gff.parse_gff(gff_file):
        info = (annotation.start, annotation.end, annotation.uid)
        try:
            annotations[annotation.seq_id].append(info)
        except KeyError:
            annotations[annotation.seq_id] = [info]

    LOG.info("Producing Dictionary Filter")
    max_size_dict = {
        seq_id: max(x[1] for x in value)
        for seq_id, value in annotations.items()
    }

    LOG.info("Finished Dictionary Filter")

    depth_data = align.SamtoolsDepth(depth_file, max_size_dict=max_size_dict, num_seqs=None)

    completed_seqs = 0

    if progress:
        pbar = tqdm(desc="Annotations", total=len(annotations))
    else:
        pbar = None

    while (completed_seqs < len(annotations)):
        seq_id = depth_data.advance_file()

        # reached the end of the depth file, exiting
        if seq_id is None:
            LOG.info("Reached the end of Depth file %r - exiting", depth_file)
            break
        
        # shouldn't be happening since we're using the max_dict_size now
        if seq_id not in annotations:
            depth_data.drop_sequence(seq_id)
            continue

        for start, end, uid in annotations[seq_id]:
            cov = depth_data.region_coverage(seq_id, start, end)
            max_cov = max(max_cov, cov)
            print(uid, "{:.5f}".format(cov), sep=separator, file=output_file)

        completed_seqs += 1
        depth_data.drop_sequence(seq_id)
        if pbar is not None:
            pbar.update(1)
    
    if pbar is not None:
        pbar.close()

    LOG.info("Found %d annotations out of %d, with a maximum coverage of %.2f", completed_seqs, len(annotations), max_cov)
