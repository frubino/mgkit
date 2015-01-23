"""
Extract information from GFF files

Sequence command
****************

Used to extract the nucleotidic sequences from GFF annotations. It requires the
*fasta* file containing the sequences referenced in the GFF `seq_id` attribute
(first column of the raw GFF).

The sequnces extract have as identifier the `uid` stored in the GFF file and by
default the sequnece is not reverse complemented if the annotation is on the
**-** strand, but this can be changed by using the `-r` option.

The sequences are wrapped at 60 characters, as per FASTA specs, but this
behavior can be disabled by specifing the `-w` option.

.. warning::

    The reference file is loaded in memory

Changes
*******

.. versionadded:: 0.1.15

"""

from __future__ import division
import sys
import argparse
import logging

import mgkit
from . import utils
from mgkit.io import gff, fasta

LOG = logging.getLogger(__name__)


def set_sequence_parser(parser):
    parser.add_argument(
        '-r',
        '--reverse',
        action='store_true',
        help='Reverse complement sequences on the - strand',
        default=False
    )
    parser.add_argument(
        '-w',
        '--no-wrap',
        action='store_true',
        help='Write the nucleotidic sequence on one line',
        default=False
    )
    parser.add_argument(
        '-f',
        '--reference',
        type=argparse.FileType('r'),
        default=None,
        help='Fasta file containing the reference sequences of the GFF file'
    )

    parser.set_defaults(func=sequence_command)


def sequence_command(options):
    if options.reference is None:
        utils.exit_script('A fasta reference file is required', 1)

    wrap = 60

    if options.no_wrap:
        wrap = None

    seqs = dict(fasta.load_fasta(options.reference))

    ann_iter = gff.parse_gff(options.input_file, gff_type=gff.from_gff)

    seq_iter = gff.extract_nuc_seqs(ann_iter, seqs, reverse=options.reverse)

    for name, seq in seq_iter:
        fasta.write_fasta_sequence(options.output_file, name, seq, wrap=wrap)


def set_common_options(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input GFF file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Output file, defaults to stdout'
    )


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Extract informations from a GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_s = subparsers.add_parser(
        'sequence',
        help='Extract the nucleotidic sequences of annotations'
    )

    set_sequence_parser(parser_s)
    set_common_options(parser_s)
    utils.add_basic_options(parser_s)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
