"""
.. versionadded:: 0.3.0

Scripts that includes some functionality to help use FASTA files with the
framework

split command
*************

Used to split a fasta file into smaller fragments

translate command
*****************

Used to translate nucleotide sequences into amino acids.

uid command
***********

Used to change a FASTA file headers to a unique ID. A table (tab separated)
with the changes made can be kept, using the *--table* option.

Changes
*******

.. versionadded:: 0.3.0

.. versionchanged:: 0.3.1
    added *translate* and *uid* command

"""

from __future__ import division
import argparse
import logging
import sys
from uuid import uuid4

import mgkit
from . import utils
from mgkit.io import fasta, open_file
from ..utils import trans_tables
from ..utils.sequence import translate_sequence

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

    parser_translate = subparsers.add_parser(
        'translate',
        help='Translate FASTA file in all 6 frames'
    )

    set_translate_parser(parser_translate)
    utils.add_basic_options(parser_translate, manual=__doc__)

    parser_uid = subparsers.add_parser(
        'uid',
        help='Changes the header to a uid (unique ID)'
    )

    set_uid_parser(parser_uid)
    utils.add_basic_options(parser_uid, manual=__doc__)

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


def set_translate_parser(parser):
    parser.add_argument(
        '-t',
        '--trans-table',
        default='universal',
        action='store',
        choices=[
            table_name.lower() for table_name in dir(trans_tables)
            if not table_name.startswith('_')
        ],
        help='translation table'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTA file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Input FASTA file, defaults to stdin'
    )
    parser.set_defaults(func=translate_command)


def load_trans_table(table_name):
    "Loads translation table "
    return getattr(trans_tables, table_name.upper())


def translate_seq(name, seq, trans_table):
    "Tranlates sequence into the 6 frames"
    header = "{0}-{1}{2}"
    for start in range(3):
        yield header.format(name, 'f', start), translate_sequence(seq, start, trans_table, False)
        yield header.format(name, 'r', start), translate_sequence(seq, start, trans_table, True)


def translate_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    trans_table = load_trans_table(options.trans_table)

    for name, seq in fasta.load_fasta(options.input_file):
        for new_header, new_seq in translate_seq(name, seq, trans_table):
            fasta.write_fasta_sequence(options.output_file, new_header, new_seq)


def set_uid_parser(parser):
    parser.add_argument(
        '-t',
        '--table',
        default=None,
        action='store',
        type=str,
        help='Filename of table of the changes (by default discards it)'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTA file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Input FASTA file, defaults to stdin'
    )
    parser.set_defaults(func=uid_command)


def uid_command(options):
    if options.table is not None:
        table_file = open_file(options.table, 'w')
        LOG.info(
            'Writing Table to file (%s)',
            getattr(options.table, 'name', repr(options.table))
        )
    else:
        table_file = None

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    for name, seq in fasta.load_fasta(options.input_file):
        uid = str(uuid4())
        if table_file is not None:
            table_file.write("{}\t{}\n".format(uid, name))

        fasta.write_fasta_sequence(options.output_file, uid, seq)


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
