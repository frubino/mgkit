"""
Commands
--------

* Interleave/deinterleave paired-end fastq files.
* Converts to FASTA
* sort 2 files to sync the headers

Changes
-------

.. versionchanged:: 0.3.4
    moved to use click, internal fastq parsing, removed *rand* command

.. versionchanged:: 0.3.1
    added stdin/stdout defaults for some commands

.. versionchanged:: 0.3.0
    added *convert* command to FASTA

"""
from __future__ import division
from builtins import zip
import logging
import click
import mgkit
from mgkit.io.fastq import choose_header_type, write_fastq_sequence, load_fastq
from mgkit.io import fasta
from . import utils

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('convert', help="""Convert FastQ file [fastq-file] to FASTA file
              [fasta-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.argument('fastq-file', type=click.File('rb'), default='-')
@click.argument('fasta-file', type=click.File('wb'), default='-')
def convert_command(verbose, fastq_file, fasta_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        "Writing FASTA file (%s)",
        getattr(fasta_file, 'name', repr(fasta_file))
    )

    for seq_id, seq, qual in load_fastq(fastq_file):
        fasta.write_fasta_sequence(fasta_file, seq_id, seq)


def report_counts(count, wcount, counter=None):
    "Logs the status"
    if (counter is None) or (counter % 100000 == 0):
        LOG.info(
            "Read %-9d sequences, Wrote %-9d (%.2f%%)",
            count,
            wcount,
            wcount / count * 100
        )


@main.command('sort', help="""Sort paired-end sequences from [mate1-input] and
              [mate2-input] into files [mate1-output] and [mate2-output]""")
@click.option('-v', '--verbose', is_flag=True)
@click.argument('mate1-input', type=click.File('rb'))
@click.argument('mate2-input', type=click.File('rb'))
@click.argument('mate1-output', type=click.File('wb'))
@click.argument('mate2-output', type=click.File('wb'))
def sort(verbose, mate1_input, mate2_input, mate1_output, mate2_output):
    "Sort two fastq files"

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing [mate1-output] to file (%s)',
        getattr(mate1_output, 'name', repr(mate1_output))
    )
    LOG.info(
        'Writing [mate2-output] to file (%s)',
        getattr(mate2_output, 'name', repr(mate2_output))
    )

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for (seq_id1, seq1, qual1), (seq_id2, seq2, qual2) in zip(load_fastq(mate1_input), load_fastq(mate2_input)):

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(seq_id1)
            if regex is None:
                simple_header = True
                LOG.info("Using a simple header structure")

        if simple_header:
            key1 = seq_id1[:-1]
            key2 = seq_id2[:-1]
        else:
            match1 = regex.search(seq_id1)
            match2 = regex.search(seq_id2)

            key1 = (
                match1.group('lane'),
                match1.group('tile'),
                match1.group('xcoord'),
                match1.group('ycoord')
            )
            key2 = (
                match2.group('lane'),
                match2.group('tile'),
                match2.group('xcoord'),
                match2.group('ycoord')
            )

        seq1 = (seq_id1, seq1, qual1)
        seq2 = (seq_id2, seq2, qual2)

        if key1 == key2:
            # if the 2
            write_fastq_sequence(mate1_output, *seq1)
            write_fastq_sequence(mate2_output, *seq2)
            wcount += 1
            report_counts(count, wcount, count)
            continue

        mate1[key1] = seq1
        mate2[key2] = seq2

        if key1 in mate2:
            write_fastq_sequence(mate1_output, *mate1[key1])
            write_fastq_sequence(mate2_output, *mate2[key1])
            del mate1[key1]
            del mate2[key1]
            wcount += 1
        if key2 in mate1:
            write_fastq_sequence(mate1_output, *mate1[key2])
            write_fastq_sequence(mate2_output, *mate2[key2])
            del mate1[key2]
            del mate2[key2]
            wcount += 1

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)


@main.command('di', help="""Deinterleave sequences from [fastq-file], into
              [mate1-file] and [mate2-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--strip', is_flag=True, default=False,
              help="Strip additional info")
@click.argument('fastq-file', type=click.File('rb'), default='-')
@click.argument('mate1-file', type=click.File('wb'))
@click.argument('mate2-file', type=click.File('wb'))
def deinterleave(verbose, strip, fastq_file, mate1_file, mate2_file):
    "Deinterleave a fastq file"

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing [mate1-file] to file (%s)',
        getattr(mate1_file, 'name', repr(mate1_file))
    )
    LOG.info(
        'Writing [mate2-file] to file (%s)',
        getattr(mate2_file, 'name', repr(mate2_file))
    )

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for seq_id, seq, qual in load_fastq(fastq_file):

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(seq_id)
            if regex is None:
                LOG.info("Using a simple header structure")
                simple_header = True

        if simple_header:
            key = seq_id[:-1]
            mate = int(seq_id[-1])
        else:
            match = regex.search(seq_id)
            key = (
                match.group('lane'),
                match.group('tile'),
                match.group('xcoord'),
                match.group('ycoord')
            )
            mate = int(match.group('mate'))

        if strip:
            sequence_name = seq_id.split('\t')[0]
        else:
            sequence_name = seq_id

        if mate == 1:
            mate1[key] = (sequence_name, seq, qual)
        else:
            mate2[key] = (sequence_name, seq, qual)

        try:
            # if sequence header in both
            seq1 = mate1[key]
            seq2 = mate2[key]
            write_fastq_sequence(mate1_file, *seq1)
            write_fastq_sequence(mate2_file, *seq2)
            wcount += 2
            del mate1[key]
            del mate2[key]
        except KeyError:
            pass

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)


@main.command('il', help="""Interleave sequences from [mate1-file] and
              [mate2-file] into [fastq-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.argument('mate1-file', type=click.File('rb'))
@click.argument('mate2-file', type=click.File('rb'))
@click.argument('fastq-file', type=click.File('wb'), default='-')
def interleave(verbose, mate1_file, mate2_file, fastq_file):
    "Interleave two fastq files"

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing interleaved [fastq-file] to file (%s)',
        getattr(fastq_file, 'name', repr(fastq_file))
    )

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for (seq_id1, seq1, qual1), (seq_id2, seq2, qual2) in zip(load_fastq(mate1_file), load_fastq(mate2_file)):

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(seq_id1)
            if regex is None:
                simple_header = True
                LOG.info("Using a simple header structure")

        if simple_header:
            key1 = seq_id1[:-1]
            key2 = seq_id2[:-1]
        else:
            match1 = regex.search(seq_id1)
            match2 = regex.search(seq_id2)

            key1 = (
                match1.group('lane'),
                match1.group('tile'),
                match1.group('xcoord'),
                match1.group('ycoord')
            )
            key2 = (
                match2.group('lane'),
                match2.group('tile'),
                match2.group('xcoord'),
                match2.group('ycoord')
            )

        seq1 = (seq_id1, seq1, qual1)
        seq2 = (seq_id2, seq2, qual2)

        if key1 == key2:
            # if the 2
            write_fastq_sequence(fastq_file, *seq1)
            write_fastq_sequence(fastq_file, *seq2)
            wcount += 1
            report_counts(count, wcount, wcount)
            continue

        mate1[key1] = seq1
        mate2[key2] = seq2

        if key1 in mate2:
            write_fastq_sequence(fastq_file, *mate1[key1])
            write_fastq_sequence(fastq_file, *mate2[key1])
            del mate1[key1]
            del mate2[key1]
            wcount += 1
        if key2 in mate1:
            write_fastq_sequence(fastq_file, *mate1[key2])
            write_fastq_sequence(fastq_file, *mate2[key2])
            del mate1[key2]
            del mate2[key2]
            wcount += 1

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)
