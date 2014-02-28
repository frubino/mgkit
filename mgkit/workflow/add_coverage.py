"""
Adds coverage information to GFF
"""
import argparse
import logging
import pysam
from . import utils
from .. io.gff import parse_gff, write_gff
from .. align import add_coverage_info
from .. import logger

LOG = logging.getLogger(__name__)


def split_sample_alg(argument):
    "Split sample/alignment option"
    try:
        sample, bam_file_name = argument.split(',', 1)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Can't get get both sample and bam file from '%s'" % argument
        )

    return (sample, bam_file_name)


def set_parser():
    "argument parser configuration"
    parser = argparse.ArgumentParser(
        description='Adds coverage information to GFF files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_file',
        action='store',
        help='GFF file to annotate',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        'output_file',
        action='store',
        help='output file with coverage information (default: GFF)',
        type=argparse.FileType('w')
    )
    parser.add_argument(
        '-a',
        '--sample-alignment',
        action='append',
        type=split_sample_alg,
        help='Number of processors to use'
    )
    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"
    options = set_parser().parse_args()

    #configs log and set log level
    logger.config_log(options.verbose)

    samples = []
    bam_files = []

    for sample, bam_file_name in options.sample_alignment:
        samples.append(sample)
        bam_files.append(pysam.Samfile(bam_file_name, 'rb'))

    if len(samples) != len(set(samples)):
        LOG.critical("There are duplicate sample names")
        return 1

    annotations = list(parse_gff(options.input_file))

    add_coverage_info(annotations, bam_files, samples)

    write_gff(annotations, options.output_file)

if __name__ == '__main__':
    main()
