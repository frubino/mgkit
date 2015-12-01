"""
Extract information from GFF files

sequence command
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

dbm command
***********

Creates a dbm DB using the `semidbm` package. The database can then be loaded
using :class:`mgkit.db.dbm.GFFDB`

mongodb command
***************

Outputs annotations in a format supported by MongoDB. More information about it
can be found in :mod:`mgkit.db.mongo`

Changes
*******

.. versionadded:: 0.1.15

.. versionadded:: 0.2.1
    *dbm* and *mongodb* commands

"""

from __future__ import division
import sys
import argparse
import logging
import functools

import mgkit
from . import utils
from mgkit.io import gff, fasta
from mgkit.db import dbm
from mgkit import taxon
from mgkit import simple_cache

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


def set_dbm_parser(parser):
    parser.add_argument(
        '-d',
        '--output-dir',
        default='gff-dbm',
        type=str,
        help='Directory for the database'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input GFF file, defaults to stdin'
    )

    parser.set_defaults(func=dbm_command)


def dbm_command(options):
    db = dbm.create_gff_dbm(
        gff.parse_gff(options.input_file),
        options.output_dir
    )
    db.close()


def set_mongodb_parser(parser):
    parser.add_argument(
        '-t',
        '--taxonomy',
        type=str,
        default=None,
        help='Taxonomy used to populate the lineage'
    )
    parser.add_argument(
        '-c',
        '--no-cache',
        action='store_false',
        default=True,
        help='No cache for the lineage function'
    )

    parser.set_defaults(func=mongodb_command)


def mongodb_command(options):

    lineage_func = None

    if options.taxonomy is not None:
        taxonomy = taxon.UniprotTaxonomy(options.taxonomy)
        lineage_func = functools.partial(
            taxon.get_lineage,
            taxonomy
        )
        if options.no_cache:
            LOG.info('Using cached calls to lineage')
            lineage_func = simple_cache.memoize(lineage_func)

    for annotation in gff.parse_gff(options.input_file):
        options.output_file.write(
            annotation.to_mongodb(lineage_func=lineage_func)
        )


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

    parser_d = subparsers.add_parser(
        'dbm',
        help='Creates a dbm database'
    )

    set_dbm_parser(parser_d)
    utils.add_basic_options(parser_d)

    parser_m = subparsers.add_parser(
        'mongodb',
        help='Extract annotations from a GFF file and makes output for MongoDB'
    )

    set_mongodb_parser(parser_m)
    set_common_options(parser_m)
    utils.add_basic_options(parser_m)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
