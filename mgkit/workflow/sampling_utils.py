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

.. versionchanged:: 0.3.4
    using *click* instead of *argparse. Now *rand_seq* can save and reload models

.. versionchanged:: 0.3.3
    added *sync*, *sample_stream* and *rand_seq* commnads

"""
from __future__ import division
from builtins import range, zip
import logging
import itertools
import uuid
import click
import numpy
import scipy.stats
import pickle
from tqdm import tqdm
import mgkit
from . import utils
from ..utils import sequence
from mgkit.io import fasta, open_file
from mgkit.io.fastq import load_fastq, write_fastq_sequence, choose_header_type

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


def infer_parameters(file_handle, fastq_bool, progress):
    LOG.info("Extrapolating model from file %s", file_handle.name)

    if fastq_bool:
        it = load_fastq(file_handle, num_qual=True)
        quals = []
    else:
        it = fasta.load_fasta(file_handle)

    if progress:
        it = tqdm(it)

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


@main.command('rand_seq', help='Generates random FastA/Q sequences')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-n', '--num-seqs', default=1000, type=click.INT, show_default=True,
              help='Number of sequences to generate')
@click.option('-gc', '--gc-content', default=.5, type=click.FLOAT, show_default=True,
              help='GC content (defaults to .5 out of 1)')
@click.option('-i', '--infer-params', default=None, type=click.File('rb'),
              help='Infer parameters GC content and Quality model from file')
@click.option('-r', '--coding-prop', default=0., type=click.FLOAT, show_default=True,
              help='Proportion of coding sequences')
@click.option('-l', '--length', default=150, type=click.INT, show_default=True,
              help='Sequence length')
@click.option('-d', '--const-model', default=False, is_flag=True,
              help='Use a model with constant qualities + noise')
@click.option('-x', '--dist-loc', default=30., type=click.FLOAT, show_default=True,
              help='Use as the starting point quality')
@click.option('-q', '--fastq', default=False, is_flag=True,
              help='The output file is a FastQ file')
@click.option('-m', '--save-model', default=None, type=click.File('wb'),
              help='Save inferred qualities model to a pickle file')
@click.option('-a', '--read-model', default=None, type=click.File('rb'),
              help='Load qualities model from a pickle file')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('output_file', type=click.File('wb'), default='-')
def rand_sequence_command(verbose, num_seqs, gc_content, infer_params,
                          coding_prop, length, const_model, dist_loc, fastq,
                          save_model, read_model, progress, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if fastq:
        if const_model:
            LOG.info("Using constant model with loc=%.1f", dist_loc)
            model = sequence.qualities_model_constant(
                length=length,
                loc=dist_loc
            )
        elif infer_params:
            length, gc_content, model = infer_parameters(
                infer_params,
                fastq,
                progress
            )
        elif read_model:
            LOG.info('Reading saved model')
            read_model = pickle.load(read_model)
            gc_content = read_model['gc_content']
            lw = read_model['lw']
            length = len(lw)
            model = (
                lw,
                getattr(scipy.stats, read_model['dist_family'])(*read_model['dist'])
            )
        else:
            LOG.info("Using decrease model with loc=%.1f", dist_loc)
            model = sequence.qualities_model_decrease(
                length=length,
                loc=dist_loc
            )

        if save_model is not None:
            LOG.info('Saving model to file (%s)',
                     getattr(save_model, 'name', repr(save_model)))
            pickle.dump(
                dict(lw=model[0], dist=model[1].args, dist_family='norm',
                     gc_content=gc_content),
                save_model
            )

    # A C T G
    prob = [
        (1 - gc_content) / 2.,
        gc_content / 2.
    ] * 2

    LOG.info(
        '%d Sequences, with a length of %d - coding proportion: %.1f',
        num_seqs,
        length,
        coding_prop
    )
    LOG.info("Probability A %.2f, C %.2f, T %.2f, G %.2f", *prob)

    num_coding = numpy.round(num_seqs * coding_prop).astype(int)
    seq_it = itertools.chain(
        sequence.random_sequences_codon(n=num_coding, length=length),
        sequence.random_sequences(
            n=num_seqs - num_coding,
            length=length,
            p=prob
        )
    )
    if fastq:
        qual_it = sequence.random_qualities(
            n=num_seqs,
            length=length,
            model=model
        )
    else:
        qual_it = itertools.repeat(num_seqs)

    if progress:
        qual_it = tqdm(qual_it, total=num_seqs)

    for seq, qual in zip(seq_it, qual_it):
        seq_id = str(uuid.uuid4())
        if fastq:
            write_fastq_sequence(output_file, seq_id, seq, qual)
        else:
            fasta.write_fasta_sequence(output_file, seq_id, seq)


def compare_header(header1, header2, header_type=None):

    if header_type is None:
        return header1[-1] == header2[-1]
    else:
        return header1.split(' ')[0] == header2.split(' ')[0]


@main.command('sync', help="""Syncs a FastQ file generated with *sample* with
              the original pair of files.""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-m', '--master-file', type=click.File('rb'), required=True,
              help='Resampled FastQ file that is out of sync with the original pair')
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def fq_sync_command(verbose, master_file, input_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    master_file = load_fastq(master_file, num_qual=False)
    master_header = next(master_file)[0]

    header_type = choose_header_type(master_header)

    written_count = 0

    for header, seq, qual in load_fastq(input_file, num_qual=False):

        if compare_header(master_header, header, header_type):
            write_fastq_sequence(output_file, header, seq, qual)
            written_count += 1
            try:
                master_header = next(master_file)[0]
            except StopIteration:
                break

    LOG.info("Wrote %d FASTQ sequences", written_count)


@main.command('sample', help='Sample a FastA/Q multiple times')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-p', '--prefix', default='sample', show_default=True,
              help='Prefix for the file name(s) in output')
@click.option('-n', '--number', type=click.INT, default=1, show_default=True,
              help='Number of samples to take')
@click.option('-r', '--prob', type=click.FLOAT, default=10**-3,
              show_default=True, help='Probability of picking a sequence')
@click.option('-x', '--max-seq', type=click.INT, default=10**5,
              show_default=True, help='Maximum number of sequences')
@click.option('-q', '--fastq', default=False, is_flag=True,
              help='The input file is a fastq file')
@click.option('-z', '--gzip', is_flag=True, default=False,
              help='gzip output files')
@click.argument('input-file', type=click.File('rb'), default='-')
def sample_command(verbose, prefix, number, prob, max_seq, fastq, gzip,
                   input_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        "Sampling %s file (%d) chunks with prefix (%s)",
        'FastQ' if fastq else 'Fasta',
        number,
        prefix
    )

    if (prob > 1) or (prob <= 0):
        utils.exit_script(
            "The probability value ({}) is outside the correct range" +
            " (0 < p <= 1)",
            1
        )

    dist = scipy.stats.binom(1, prob)

    LOG.info(
        "Probability of picking a sequence (%.5f), max number of seqs %d",
        prob,
        max_seq
    )
    name_mask = "%s-{0:05}.%s" % (prefix, 'fq' if fastq else 'fa')

    if gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    output_files = [
        dict(
            h=open_file(name_mask.format(i), 'wb'),
            c=0
        )
        for i in range(number)
    ]

    load_func = load_fastq if fastq else fasta.load_fasta
    write_func = write_fastq_sequence if fastq else fasta.write_fasta_sequence

    for seq in load_func(input_file):
        # reached the maximum number of sequences for all samples
        if all(x['c'] == max_seq for x in output_files):
            break

        for output in output_files:
            if output['c'] == max_seq:
                continue

            if dist.rvs():
                write_func(output['h'], *seq)
                output['c'] += 1


@main.command('sample_stream', help="""Samples a FastA/Q one time, alternative
              to sample if multiple sampling is not needed""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-r', '--prob', type=click.FLOAT, default=10**-3,
              help='Probability of picking a sequence')
@click.option('-x', '--max-seq', type=click.INT, default=10**5,
              help='Maximum number of sequences')
@click.option('-q', '--fastq', default=False, is_flag=True,
              help='The input file is a fastq file')
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def stream_sample_command(verbose, prob, max_seq, fastq, input_file,
                          output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if (prob > 1) or (prob <= 0):
        utils.exit_script(
            "The probability value ({}) is outside the correct range" +
            " (0 < p <= 1)",
            1
        )

    dist = scipy.stats.binom(1, prob)

    LOG.info(
        "Probability of picking a sequence (%.5f), max number of seqs %d",
        prob,
        max_seq
    )

    load_func = load_fastq if fastq else fasta.load_fasta
    write_func = write_fastq_sequence if fastq else fasta.write_fasta_sequence

    count = 0

    for index, seq in enumerate(load_func(input_file)):
        # reached the maximum number of sequences for all samples
        if count >= max_seq:
            LOG.info('Read the first %d sequences', index + 1)
            break

        if dist.rvs():
            write_func(output_file, *seq)
            count += 1
