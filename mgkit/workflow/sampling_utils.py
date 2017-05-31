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

    return parser


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

    if (options.prob == 1) or (options.prob <= 0):
        utils.exit_script(
            "The probability value ({}) is outside the correct range" +
            " (0 < p < 1)",
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
