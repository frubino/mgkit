"""
Interleave/deinterleave paired-end fastq files.
"""
from __future__ import division

import random
import argparse
import logging
import mgkit
import HTSeq
from itertools import izip
from mgkit.io.fastq import choose_header_type
from mgkit.io.fastq import write_fastq_sequence
from . import utils

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Interleave/Deinterleave fastq files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    utils.add_basic_options(parser)
    #Deinterleave
    subparsers = parser.add_subparsers()
    parser_d = subparsers.add_parser('di', help='Deinterleave sequences')
    parser_d.add_argument(
        'input_file',
        type=argparse.FileType('r'),
        action='store',
        help="File with interleaved sequences"
    )
    parser_d.add_argument(
        'mate1_file',
        type=argparse.FileType('w'),
        action='store',
        help="File with mate 1 reads"
    )
    parser_d.add_argument(
        'mate2_file',
        type=argparse.FileType('w'),
        action='store',
        help="File with mate 2 reads"
    )
    parser_d.add_argument(
        '-s',
        '--strip',
        action='store_true',
        default=False,
        help="Strip additional info"
    )
    parser_d.set_defaults(func=deinterleave)
    #Interleave
    parser_i = subparsers.add_parser('il', help='Interleave sequences')
    parser_i.add_argument(
        'mate1_file',
        type=argparse.FileType('r'),
        action='store',
        help="File with mate 1 reads"
    )
    parser_i.add_argument(
        'mate2_file',
        type=argparse.FileType('r'),
        action='store',
        help="File with mate 1 reads"
    )
    parser_i.add_argument(
        'output_file',
        default='output-il.fq',
        type=argparse.FileType('w'),
        action='store',
        help='Output file name'
    )
    parser_i.set_defaults(func=interleave)
    #Sort
    parser_s = subparsers.add_parser('sort', help='Sort paired-end sequences')
    parser_s.add_argument(
        'mate1_input',
        type=argparse.FileType('r'),
        action='store',
        help="File with mate 1 reads"
    )
    parser_s.add_argument(
        'mate2_input',
        type=argparse.FileType('r'),
        action='store',
        help="File with mate 1 reads"
    )
    parser_s.add_argument(
        'mate1_output',
        type=argparse.FileType('w'),
        action='store',
        help="Output file with mate 1 reads"
    )
    parser_s.add_argument(
        'mate2_output',
        type=argparse.FileType('w'),
        action='store',
        help="Output file with mate 2 reads"
    )
    parser_s.set_defaults(func=sort)
    #Randomise
    parser_r = subparsers.add_parser('rand', help='Randomise sequences')
    parser_r.add_argument(
        'input_file',
        type=argparse.FileType('r'),
        action='store',
        help="File with sequences"
    )
    parser_r.add_argument(
        'output_file',
        type=argparse.FileType('w'),
        action='store',
        help="Output file"
    )
    parser_r.add_argument(
        '-b',
        '--buffer-size',
        type=int,
        action='store',
        default=None,
        help="Maximum number of sequences to keep in memory (default: all)"
    )
    parser_r.set_defaults(func=randomise)

    return parser


def report_counts(count, wcount, counter=None):
    "Logs the status"
    if (counter is None) or (counter % 100000 == 0):
        LOG.info(
            "Read %-7d sequences, Wrote %-7d (%.2f%%)",
            count,
            wcount,
            wcount / count * 100
        )


def randomise(options):
    "Used for testing"
    LOG.info(
        "Reading from file (%s)",
        options.input_file.name
    )

    file_handle = HTSeq.FastqReader(options.input_file)
    random.seed(random.randint(0, 1000000))

    count = 0
    wcount = 0

    seqs = []

    for sequence in file_handle:
        count += 1
        report_counts(count, wcount, count)
        seqs.append((sequence.name, sequence.seq, sequence.qualstr))

        if (options.buffer_size is None) or (len(seqs) < options.buffer_size):
            continue

        index = random.randint(0, len(seqs) - 1)

        write_fastq_sequence(options.output_file, *seqs[index])

        del seqs[index]

        wcount += 1

    #Shuffle the list
    random.shuffle(seqs)

    #Flush the rest of the sequences to disk
    for sequence in seqs:
        write_fastq_sequence(options.output_file, *sequence)
        wcount += 1
        report_counts(count, wcount, wcount)

    report_counts(count, wcount, None)


def sort(options):
    "Sort two fastq files"
    LOG.info(
        "Reading from mate1 (%s) and mate2 (%s)",
        options.mate1_input.name,
        options.mate2_input.name
    )

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for sequence1, sequence2 in izip(HTSeq.FastqReader(options.mate1_input),
                                     HTSeq.FastqReader(options.mate2_input)):

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(sequence1.name)
            if regex is None:
                simple_header = True
                LOG.info("Using a simple header structure")

        if simple_header:
            key1 = sequence1.name[:-1]
            key2 = sequence2.name[:-1]
        else:
            match1 = regex.search(sequence1.name)
            match2 = regex.search(sequence2.name)

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

        seq1 = (sequence1.name, sequence1.seq, sequence1.qualstr)
        seq2 = (sequence2.name, sequence2.seq, sequence2.qualstr)

        if key1 == key2:
            #if the 2
            write_fastq_sequence(options.mate1_output, *seq1)
            write_fastq_sequence(options.mate2_output, *seq2)
            wcount += 1
            report_counts(count, wcount, count)
            continue

        mate1[key1] = seq1
        mate2[key2] = seq2

        if key1 in mate2:
            write_fastq_sequence(options.mate1_output, *mate1[key1])
            write_fastq_sequence(options.mate2_output, *mate2[key1])
            del mate1[key1]
            del mate2[key1]
            wcount += 1
        if key2 in mate1:
            write_fastq_sequence(options.mate1_output, *mate1[key2])
            write_fastq_sequence(options.mate2_output, *mate2[key2])
            del mate1[key2]
            del mate2[key2]
            wcount += 1

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)


def deinterleave(options):
    "Deinterleave a fastq file"

    LOG.info("Reading from file %s", options.input_file.name)

    fastq_in = HTSeq.FastqReader(options.input_file)

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for sequence in fastq_in:

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(sequence.name)
            if regex is None:
                LOG.info("Using a simple header structure")
                simple_header = True

        if simple_header:
            key = sequence.name[:-1]
            mate = int(sequence.name[-1])
        else:
            match = regex.search(sequence.name)
            key = (
                match.group('lane'),
                match.group('tile'),
                match.group('xcoord'),
                match.group('ycoord')
            )
            mate = int(match.group('mate'))

        if options.strip:
            sequence_name = sequence.name.split('\t')[0]
        else:
            sequence_name = sequence.name

        if mate == 1:
            mate1[key] = (sequence_name, sequence.seq, sequence.qualstr)
        else:
            mate2[key] = (sequence_name, sequence.seq, sequence.qualstr)

        try:
            #if sequence header in both
            seq1 = mate1[key]
            seq2 = mate2[key]
            write_fastq_sequence(options.mate1_file, *seq1)
            write_fastq_sequence(options.mate2_file, *seq2)
            wcount += 2
            del mate1[key]
            del mate2[key]
        except KeyError:
            pass

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)


def interleave(options):
    "Interleave two fastq files"

    LOG.info(
        "Reading from mate1 (%s) and mate2 (%s)",
        options.mate1_file.name,
        options.mate2_file.name
    )

    regex = None
    simple_header = False

    mate1 = {}
    mate2 = {}

    count = 0
    wcount = 0

    for sequence1, sequence2 in izip(HTSeq.FastqReader(options.mate1_file),
                                     HTSeq.FastqReader(options.mate2_file)):

        count += 1

        if (regex is None) and (not simple_header):
            regex = choose_header_type(sequence1.name)
            if regex is None:
                simple_header = True
                LOG.info("Using a simple header structure")

        if simple_header:
            key1 = sequence1.name[:-1]
            key2 = sequence2.name[:-1]
        else:
            match1 = regex.search(sequence1.name)
            match2 = regex.search(sequence2.name)

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

        seq1 = (sequence1.name, sequence1.seq, sequence1.qualstr)
        seq2 = (sequence2.name, sequence2.seq, sequence2.qualstr)

        if key1 == key2:
            #if the 2
            write_fastq_sequence(options.output_file, *seq1)
            write_fastq_sequence(options.output_file, *seq2)
            wcount += 1
            report_counts(count, wcount, wcount)
            continue

        mate1[key1] = seq1
        mate2[key2] = seq2

        if key1 in mate2:
            write_fastq_sequence(options.output_file, *mate1[key1])
            write_fastq_sequence(options.output_file, *mate2[key1])
            del mate1[key1]
            del mate2[key1]
            wcount += 1
        if key2 in mate1:
            write_fastq_sequence(options.output_file, *mate1[key2])
            write_fastq_sequence(options.output_file, *mate2[key2])
            del mate1[key2]
            del mate2[key2]
            wcount += 1

        report_counts(count, wcount, count)

    report_counts(count, wcount, None)


def main():
    "Main function"
    options = set_parser().parse_args()

    #configs log and set log level
    mgkit.logger.config_log(options.verbose)
    options.func(options)

if __name__ == '__main__':
    main()
