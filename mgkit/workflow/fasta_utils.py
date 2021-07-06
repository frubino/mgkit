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

filter
******

Used to filter a FASTA file by length and also for sequence/header if a pattern is contained.
A list of headers to keep can be passed using the `-f` option.

info
****

Gets information about a FASTA file, prints seq_id (trimmed at first space), length
and hash (default sha1) and optionally the sequence, GC content and in GFF format if wanted.

rename
******

Renames the headers of a FASTA file, appending a random suffix and an optional prefix

Changes
*******

.. versionadded:: 0.3.0

.. versionchanged:: 0.3.1
    added *translate* and *uid* command

.. versionchanged:: 0.3.4
    ported to *click*

.. versionchanged:: 0.5.5
    added option `-1` to output only the forward/frame0 and `-w` to avoid wrap at 60 chars to the *translate* command

.. versionchanged:: 0.5.7
    added `filter` and `info` commands for simple fasta file filtering and info

"""

from __future__ import division
from builtins import range
import logging
from uuid import uuid4
import click
import hashlib
import pathlib
import string
import random
from tqdm import tqdm
import mgkit
from . import utils
from mgkit.io import fasta, gff
from ..utils import trans_tables
from ..utils.sequence import translate_sequence, sequence_gc_content

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
@click.option('-1', '--one-seq', default=False, is_flag=True, show_default=True,
              help='Only translate the sequence, instead of all 6 frames')
@click.option('-w', '--no-wrap', default=False, is_flag=True, show_default=True,
              help='Make a sequence use only 1 line (2 including header)')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def translate_command(verbose, trans_table, one_seq, no_wrap, progress, fasta_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    if one_seq:
        LOG.info("Assuming the sequences are in the correct frame")
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    trans_table = load_trans_table(trans_table)

    iterator = fasta.load_fasta(fasta_file)

    if progress:
        iterator = tqdm(iterator)
    
    if no_wrap:
        wrap = None
    else:
        wrap = 60

    for name, seq in iterator:
        if one_seq:
            new_seq = translate_sequence(seq, 0, trans_table, False)
            fasta.write_fasta_sequence(output_file, name, new_seq, wrap=wrap)
        else:
            for new_header, new_seq in translate_seq(name, seq, trans_table):
                fasta.write_fasta_sequence(output_file, new_header, new_seq, wrap=wrap)


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


@main.command('filter', help="""Filters a FASTA file [file-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('--len-gt', default=None, type=click.IntRange(min=1), help='Keeps sequences whose length is greater than')
@click.option('--len-lt', default=None, type=click.IntRange(min=1), help='Keeps sequences whose length is less than')
@click.option('--header-contains', default=None, type=click.STRING, help='Keeps sequences whose header contains the string')
@click.option('--seq-pattern', default=None, type=click.STRING, help='Keeps sequences that contains the string')
@click.option('-f', '--header-file', default=None, type=click.File('r'), 
                help='Keep only sequences contained in file list')
@click.option('-w', '--wrap', default=False, is_flag=True, help='Wraps the output sequences to 60 characters')
@click.option('-s', '--trim-tail', default=False, is_flag=True,
                help='Removes header information after first space')
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def filter_command(verbose, len_gt, len_lt, header_contains, seq_pattern, header_file, wrap, 
                    trim_tail, fasta_file, output_file):
    """
    .. versionadded:: 0.5.7

    Filters a fasta file
    """
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    
    if wrap:
        wrap = 60
    else:
        wrap = None
    
    if header_file is not None:
        LOG.info('Reading header list from file %r', header_file)
        header_file = set(
            line.strip()
            for line in header_file
        )

    count = 0

    for name, seq in fasta.load_fasta(fasta_file):
        seq_len = len(seq)
        if len_gt is not None:
            if not (seq_len > len_gt):
                continue
        if len_lt is not None:
            if not (seq_len < len_lt):
                continue
        
        if trim_tail:
            name = name.split(' ')[0]

        if header_file is not None:
            if name not in header_file:
                continue

        if header_contains is not None:
            if header_contains not in name:
                continue
        if seq_pattern is not None:
            if seq_pattern not in seq:
                continue
        
        fasta.write_fasta_sequence(output_file, name, seq, wrap=wrap)
        count += 1

    LOG.info('Kept %d sequence(s)', count)


@main.command('info', help="""Gets information of FASTA file [file-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-h', '--header', is_flag=True, default=False,
    help="Prints header")
@click.option('-s', '--include-seq', is_flag=True, default=False,
    help="Includes the sequence")
@click.option('-r', '--no-rename', is_flag=True, default=False,
    help="Do not split sequence name at first space")
@click.option('-a', '--hash-type', type=click.Choice(['sha1', 'md5', 'sha256']),
    default='sha1', show_default=True)
@click.option('-g', '--out-gff', is_flag=True, default=False, show_default=True,
    help='Outputs a GFF file')
@click.option('-gc', '--gc-content', is_flag=True, default=False, show_default=True,
    help='Includes the GC Content')
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('w'), default='-')
def info_command(verbose, header, include_seq, no_rename, hash_type, out_gff, gc_content, fasta_file, output_file):
    """
    .. versionadded:: 0.5.7

    Makes a tabular format file with hash, length and optionally the sequence. Alternatively, a GFF
    can be the output.
    """
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    hash_func = getattr(hashlib, hash_type)
    if no_rename:
        LOG.info("Sequence names will not be split at the first space")
        load_func = fasta.load_fasta
    else:
        load_func = fasta.load_fasta_rename

    if header and (not out_gff):
        columns = ["seq_id", "length", "hash"]
        if include_seq:
            columns.append('sequence')
        if gc_content:
            columns.append('gc_content')
        print(*columns, sep='\t', file=output_file)
    
    for name, seq in load_func(fasta_file):
        seq_hash = hash_func(seq.encode('ascii')).hexdigest()
        if gc_content:
            gc_value = sequence_gc_content(seq)

        if out_gff:
            annotation = gff.from_sequence(name, seq, seq_hash=seq_hash)
            if gc_content:
                annotation.set_attr("gc_content", gc_value)
            output_file.write(annotation.to_gff())
        else:
            values = [name, len(seq), hash_func(seq.encode('ascii')).hexdigest()]
            if include_seq:
                values.append(seq)
            if gc_content:
                values.append(gc_value)
            print(*values, sep='\t', file=output_file)


@main.command('rename', help="""Rename Sequence headers of FASTA file [file-file]
                Adds 2 possible elements to the sequence header, separated by a
                character 1) a suffix (random string of characters) and 2) a prefix
                (optional).

                The character used as separator should be a '|' (default), '#' or other
                character that is not truncated in other software (space is).

                In fact, this script will truncate the header at the first space
                """)
@click.option('-v', '--verbose', is_flag=True)
@click.option('-p', '--prefix', default=None, type=click.STRING,
              help="Adds a prefix to the header")
@click.option('-f', '--file-name', is_flag=True, default=False,
              help="Adds filename as prefix (Useful for adding the file name")
@click.option('-s', '--separator', default='|', type=click.STRING,
              help="Separator for the elements of the new header")
@click.option('-l', '--suffix-len', default=5, type=click.IntRange(0, 20),
              help="Number of random characters to use")
@click.argument('fasta-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def rename_command(verbose, prefix, file_name, separator, suffix_len, fasta_file, output_file):
    """
    .. versionadded:: 0.5.7

    Renames headers of a fasta file
    """
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    random_pop = string.ascii_letters + string.digits
    if prefix is not None:
        if ' ' in prefix:
            utils.exit_script("Prefix contains spaces: `{}`".format(prefix), 1)
    if file_name:
        file_name = pathlib.Path(fasta_file.name).name
        if file_name == '<stdin>':
            utils.exit_script("Standard input used, cannot use `-f`", 2)

    base_elements = []
    if prefix is not None:
        base_elements.append(prefix)
    if file_name:
        base_elements.append(file_name)

    for name, seq in fasta.load_fasta_rename(fasta_file):
        elements = base_elements + [name]
        if suffix_len > 0:
            elements.append(
                ''.join(random.sample(random_pop, k=suffix_len))
            )
        fasta.write_fasta_sequence(output_file, separator.join(elements), seq)
