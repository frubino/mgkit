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

*sample_stream* command
-----------------------

It works in the same way as *sample*, however the file is sampled only once and
the output is the stdout by default. This can be convenient if streams are a
preferred way to sample the file.

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
    it's unusual to find it in SRA or other repositories as well.

*rand_seq* command
------------------

Generate random FastA/Q sequences, allowing the specification of GC content and
number of sequences being coding or random. If the output format chosen is
FastQ, qualities are generated using a decreasing model with added noise. A
constant model can be specified instead with a switch. Parameters such GC,
length and the type of model can be infered by passing a FastA/Q file, with
the quality model fit using a LOWESS (using :func:`mgkit.utils.sequence.extrapolate_model`).
The noise in that case is model as the a normal distribution fitted from the
qualities along the sequence deviating from the fitted LOWSS and scaled back by
half to avoid too drastic changes in the qualities. Also the qualities are
clipped at 40 to avoid compatibility problems with FastQ readers. If inferred,
the model can be saved (as a pickle file) and loaded back for analysis

Changes
-------

.. versionchanged:: 0.3.3
    added *sync*, *sample_stream* and *rand_seq* commnads

"""
from __future__ import division
import argparse
import logging
import itertools
import uuid
import numpy
import scipy.stats
import pickle

import mgkit
from . import utils
from ..utils import sequence
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
        help='Sample a FastA/Q multiple times'
    )

    set_sample_parser(parser_sample)
    utils.add_basic_options(parser_sample, manual=__doc__)

    parser_fq_sync = subparsers.add_parser(
        'sync',
        help='Syncs a FastQ'
    )

    set_fq_sync_parser(parser_fq_sync)
    utils.add_basic_options(parser_fq_sync, manual=__doc__)

    stream_sample_parser = subparsers.add_parser(
        'sample_stream',
        help='Sample a FastA/Q one time'
    )

    set_stream_sample_parser(stream_sample_parser)
    utils.add_basic_options(stream_sample_parser, manual=__doc__)

    random_seq_parser = subparsers.add_parser(
        'rand_seq',
        help='Generates random FastA/Q sequences'
    )

    set_random_seq_parser(random_seq_parser)
    utils.add_basic_options(random_seq_parser, manual=__doc__)

    return parser


def set_random_seq_parser(parser):
    parser.add_argument(
        '-n',
        '--num-seqs',
        default=1000,
        type=int,
        help='Number of sequences to generate'
    )
    parser.add_argument(
        '-gc',
        '--gc-content',
        default=.5,
        type=float,
        help='GC content (defaults to .5 out of 1)'
    )
    parser.add_argument(
        '-i',
        '--infer-params',
        default=None,
        type=argparse.FileType('r'),
        help='Infer parameters GC content and Quality model from file'
    )
    parser.add_argument(
        '-r',
        '--coding-prop',
        default=0.,
        type=float,
        help='Proportion of coding sequences'
    )
    parser.add_argument(
        '-l',
        '--length',
        default=150,
        type=int,
        help='Sequence length'
    )
    parser.add_argument(
        '-d',
        '--const-model',
        default=False,
        action='store_true',
        help='Use a model with constant qualities + noise'
    )
    parser.add_argument(
        '-x',
        '--dist-loc',
        default=30.,
        type=float,
        help='Use as the starting point quality'
    )
    parser.add_argument(
        '-q',
        '--fastq',
        default=False,
        action='store_true',
        help='The output file is a FastQ file'
    )
    parser.add_argument(
        '-m',
        '--save-model',
        default=None,
        type=str,
        help='Save inferred qualities model to a pickle file'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default='-',
        help='Output FastA/Q file, defaults to stdout'
    )
    parser.set_defaults(func=rand_sequence_command)


def infer_parameters(file_handle, fastq_bool):
    LOG.info('Inferring parameters from file')

    if fastq_bool:
        it = fastq.load_fastq(file_handle, num_qual=True)
        quals = []
    else:
        it = fasta.load_fasta(file_handle)

    gc_content = []

    length = 0

    for record in it:
        length = max(length, len(record[1]))
        gc_content.append(
            sequence.sequence_gc_content(record[1])
        )
        if fastq_bool:
            quals.append(record[2])

    if fastq_bool:
        model = sequence.extrapolate_model(quals)
    else:
        model = None

    gc_content = numpy.mean(gc_content)

    return length, gc_content, model


def rand_sequence_command(options):

    if options.infer_params:
        length, gc_content, model = infer_parameters(
            options.infer_params,
            options.fastq
        )
        if options.save_model is not None:
            pickle.dump(dict(lw=model[0], dist=model[1].args, dist_family='norm'), open_file(options.save_model, 'w'))
    else:
        model = None
        length = options.length
        gc_content = options.gc_content

    # A C T G
    prob = [
        (1 - gc_content) / 2.,
        gc_content / 2.
    ] * 2

    LOG.info(
        '%d Sequences, with a length of %d - coding proportion: %.1f',
        options.num_seqs,
        length,
        options.coding_prop
    )
    LOG.info("Probability A %.2f, C %.2f, T %.2f, G %.2f", *prob)

    if (not options.infer_params) and options.fastq:
        if options.const_model:
            LOG.info("Using constant model with loc=%.1f", options.dist_loc)
            model = sequence.qualities_model_constant(
                length=length,
                loc=options.dist_loc
            )
        else:
            LOG.info("Using decrease model with loc=%.1f", options.dist_loc)
            model = sequence.qualities_model_decrease(
                length=length,
                loc=options.dist_loc
            )

    num_coding = numpy.round(options.num_seqs * options.coding_prop).astype(int)
    seq_it = itertools.chain(
        sequence.random_sequences_codon(n=num_coding, length=length),
        sequence.random_sequences(
            n=options.num_seqs - num_coding,
            length=length,
            p=prob
        )
    )
    if options.fastq:
        qual_it = sequence.random_qualities(
            n=options.num_seqs,
            length=length,
            model=model
        )
    else:
        qual_it = itertools.repeat(options.num_seqs)

    for seq, qual in itertools.izip(seq_it, qual_it):
        seq_id = str(uuid.uuid4())
        if options.fastq:
            fastq.write_fastq_sequence(options.output_file, seq_id, seq, qual)
        else:
            fasta.write_fasta_sequence(options.output_file, seq_id, seq)


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
        '-z',
        '--gzip',
        action='store_true',
        default=False,
        help='gzip output files'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input FASTA/FASTQ file, defaults to stdin'
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


def set_stream_sample_parser(parser):
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
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default='-',
        help='Output FASTA/FASTQ file, defaults to stdout'
    )
    parser.set_defaults(func=stream_sample_command)


def stream_sample_command(options):

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

    load_func = fastq.load_fastq if options.fastq else fasta.load_fasta
    write_func = fastq.write_fastq_sequence if options.fastq else fasta.write_fasta_sequence

    count = 0

    for index, seq in enumerate(load_func(options.input_file)):
        # reached the maximum number of sequences for all samples
        if count >= options.max_seq:
            LOG.info('Read the first %d sequences', index + 1)
            break

        if dist.rvs():
            write_func(options.output_file, *seq)
            count += 1


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
