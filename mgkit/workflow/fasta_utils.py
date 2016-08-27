"""
.. versionadded:: 0.3.0

Scripts that includes some functionality to help use FASTA files with the
framework

split command
*************

Used to split a fasta file into smaller fragments

Changes
*******

.. versionadded:: 0.3.0

"""

from __future__ import division
import argparse
import logging

import mgkit
from . import utils
from mgkit.io import fasta

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='FASTA utilities',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    utils.add_basic_options(parser, manual=__doc__)

    subparsers = parser.add_subparsers()

    parser_split = subparsers.add_parser(
        'split',
        help='Splits a FASTA file in a number of fragments'
    )

    set_split_parser(parser_split)
    utils.add_basic_options(parser_split, manual=__doc__)

    return parser


def set_split_parser(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTA file, defaults to stdin'
    )
    parser.add_argument(
        '-p',
        '--prefix',
        type=str,
        default='split',
        help='Prefix for the file name in output'
    )
    parser.add_argument(
        '-n',
        '--number',
        type=int,
        default=10,
        help='Number of chunks into which split the FASTA file'
    )
    parser.add_argument(
        '-z',
        '--gzip',
        action='store_true',
        default=False,
        help='gzip output files'
    )
    parser.set_defaults(func=split_command)


def split_command(options):
    LOG.info(
        "Splitting FASTA into %d chunks with prefix %s",
        options.number,
        options.prefix
    )

    name_mask = "%s-{0:05}.fa" % options.prefix
    if options.gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    fasta.split_fasta_file(
        options.input_file,
        name_mask,
        num_files=options.number
    )


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
