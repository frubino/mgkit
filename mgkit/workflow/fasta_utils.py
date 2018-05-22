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

.. versionchanged:: 0.3.4
    ported to *click*

"""

from __future__ import division
from builtins import range
import logging
import sys
from uuid import uuid4

import click
import mgkit
from . import utils
from mgkit.io import fasta, open_file
from ..utils import trans_tables
from ..utils.sequence import translate_sequence

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('split', help="""Splits a FASTA file [fasta-file] in a number of
              fragments""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-p', '--prefix', default='split', show_default=True,
              help='Prefix for the file name in output')
@click.option('-n', '--number', type=click.INT, default=10, show_default=True,
              help='Number of chunks into which split the FASTA file')
@click.option('-z', '--gzip', is_flag=True, default=False,
              help='gzip output files')
@click.argument('fasta-file', type=click.File('rb'), default='-')
def split_command(verbose, prefix, number, gzip, fasta_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        "Splitting FASTA into %d chunks with prefix (%s)",
        number,
        prefix
    )

    name_mask = "%s-{0:05}.fa" % prefix
    if gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    fasta.split_fasta_file(
        fasta_file,
        name_mask,
        num_files=number
    )


def load_trans_table(table_name):
    "Loads translation table "
    return getattr(trans_tables, table_name.upper())


def translate_seq(name, seq, trans_table):
    "Tranlates sequence into the 6 frames"
    header = "{0}-{1}{2}"
    for start in range(3):
        yield header.format(name, 'f', start), translate_sequence(seq, start, trans_table, False)
        yield header.format(name, 'r', start), translate_sequence(seq, start, trans_table, True)


@main.command('translate', help="""Translate FASTA file [fasta-file] in all 6
              frames to [output-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--trans-table', default='universal', show_default=True,
              type=click.Choice([table_name.lower() for table_name in dir(trans_tables) if not table_name.startswith('_')]),
              help='translation table')
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def translate_command(verbose, trans_table, fasta_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    trans_table = load_trans_table(trans_table)

    for name, seq in fasta.load_fasta(fasta_file):
        for new_header, new_seq in translate_seq(name, seq, trans_table):
            fasta.write_fasta_sequence(output_file, new_header, new_seq)


@main.command('uid', help="""Changes each header of a FASTA file [file-file] to
              a uid (unique ID)""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--table', default=None, type=click.File('wb'),
              help='Filename of a table to record the changes (by default discards it)')
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def uid_command(verbose, table, fasta_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    if table is not None:
        LOG.info(
            'Writing Table to file (%s)',
            getattr(table, 'name', repr(table))
        )

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    for name, seq in fasta.load_fasta(fasta_file):
        uid = str(uuid4())
        if table is not None:
            table.write("{}\t{}\n".format(uid, name).encode('ascii'))

        fasta.write_fasta_sequence(output_file, uid, seq)
