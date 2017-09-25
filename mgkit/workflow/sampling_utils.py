"""
.. versionadded:: 0.3.1

Resampling Utilities
====================

*sample* command
----------------

This command samples from a Fasta or FastQ file, based on a probability defined
by the user (0.001 or 1 / 1000 by default, *-r* parameter), for a maximum number
of sequences (100,000 by default, *-x* parameter). By default 1 sample is
extracted, but as many as desired can be taken, by using the *-n* parameter.

The sequence file in input can be either be passed to the standard input or as
last parameter on the command line. By defult a Fasta is expected, unless the
*-q* parameter is passed.

The *-p* parameter specifies the prefix to be used, and if the output files can
be gzipped using the *-z* parameter.

*sync* command
--------------

Used to keep in sync forward and reverse read files in paired-end FASTQ.
The scenario is that the *sample* command was used to resample a FASTQ file,
usually the forward, but we need the reverse as well. In this case, the resampled
file, called *master* is passed to the *-m* option and the input file is
the file that is to be synced (reverse). The input file is scanned until the same header is
found in the master file and when that happens, the sequence is written. The
next sequence is then read from the master file and the process is repeated until all
sequence in the master file are found in the input file. This implies having
the 2 files sorted in the same way, which is what the *sample* command does.

.. note::

    the old casava format is not supported by this command at the moment, as
    it's unusual to find it in SRA or other repository as well.

Changes
-------

.. versionchanged:: 0.3.3
    added *sync* commnad

"""
from __future__ import division
import argparse
import logging
import scipy.stats

import mgkit
from . import utils
from mgkit.io import fasta, fastq, open_file

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Sampling utilities',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    utils.add_basic_options(parser, manual=__doc__)

    subparsers = parser.add_subparsers()

    parser_sample = subparsers.add_parser(
        'sample',
        help='Sample a Fasta/FastQ'
    )

    set_sample_parser(parser_sample)
    utils.add_basic_options(parser_sample, manual=__doc__)

    parser_fq_sync = subparsers.add_parser(
        'sync',
        help='Syncs a FASTQ'
    )

    set_fq_sync_parser(parser_fq_sync)
    utils.add_basic_options(parser_fq_sync, manual=__doc__)

    return parser


def set_fq_sync_parser(parser):
    parser.add_argument(
        '-m',
        '--master-file',
        type=argparse.FileType('r'),
        required=True,
        help=''
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTQ file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default='-',
        help='Input FASTQ file, defaults to stdout'
    )
    parser.set_defaults(func=fq_sync_command)


def compare_header(header1, header2, header_type=None):

    if header_type is None:
        return header1[-1] == header2[-1]
    else:
        return header1.split(' ')[0] == header2.split(' ')[0]

def fq_sync_command(options):

    master_file = fastq.load_fastq(options.master_file, num_qual=False)
    master_header = next(master_file)[0]

    header_type = fastq.choose_header_type(master_header)

    written_count = 0

    for header, seq, qual in fastq.load_fastq(options.input_file, num_qual=False):

        if compare_header(master_header, header, header_type):
            fastq.write_fastq_sequence(
                options.output_file, header, seq, qual
            )
            written_count += 1
            try:
                master_header = next(master_file)[0]
            except StopIteration:
                break

    LOG.info("Wrote %d FASTQ sequences", written_count)

def set_sample_parser(parser):
    parser.add_argument(
        '-p',
        '--prefix',
        type=str,
        default='sample',
        help='Prefix for the file name(s) in output'
    )
    parser.add_argument(
        '-n',
        '--number',
        type=int,
        default=1,
        help='Number of samples to take'
    )
    parser.add_argument(
        '-r',
        '--prob',
        type=float,
        default=10**-3,
        help='Probability of picking a sequence'
    )
    parser.add_argument(
        '-x',
        '--max-seq',
        type=int,
        default=10**5,
        help='Maximum number of sequences'
    )
    parser.add_argument(
        '-q',
        '--fastq',
        default=False,
        action='store_true',
        help='The input file is a fastq file'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTA file, defaults to stdin'
    )
    parser.add_argument(
        '-z',
        '--gzip',
        action='store_true',
        default=False,
        help='gzip output files'
    )
    parser.set_defaults(func=sample_command)


def sample_command(options):
    LOG.info(
        "Sampling %s file (%d) chunks with prefix (%s)",
        'FastQ' if options.fastq else 'Fasta',
        options.number,
        options.prefix
    )

    if (options.prob > 1) or (options.prob <= 0):
        utils.exit_script(
            "The probability value ({}) is outside the correct range" +
            " (0 < p <= 1)",
            1
        )

    dist = scipy.stats.binom(1, options.prob)

    LOG.info(
        "Probability of picking a sequence (%.5f), max number of seqs %d",
        options.prob,
        options.max_seq
    )
    name_mask = "%s-{0:05}.%s" % (options.prefix, 'fq' if options.fastq else 'fa')

    if options.gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    output_files = [
        dict(
            h=open_file(name_mask.format(i), 'w'),
            c=0
        )
        for i in xrange(options.number)
    ]

    load_func = fastq.load_fastq if options.fastq else fasta.load_fasta
    write_func = fastq.write_fastq_sequence if options.fastq else fasta.write_fasta_sequence

    for seq in load_func(options.input_file):
        # reached the maximum number of sequences for all samples
        if all(map(lambda x: x['c'] == options.max_seq, output_files)):
            break

        for output in output_files:
            if output['c'] == options.max_seq:
                continue

            if dist.rvs():
                write_func(output['h'], *seq)
                output['c'] += 1


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
