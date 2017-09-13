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

gtf command
***********

Outputs annotations in the GTF format

split command
*************

Splits a GFF file into smaller chunks, ensuring that all of a sequence
annotations are in the same file.

cov command
***********

Calculate annotation coverage for each contig in a GFF file. The command can be
run as strand specific (not by default) and requires the reference file to
which the annotation refer to. The output file is a tab separated one, with the
first column being the sequence name, the second is the strand (+, -, or NA if
not strand specific) and the third is the percentage of the sequence covered by
annotations.

.. warning::

    The GFF file is assumed to be sorted, by sequence or sequence-strand if
    wanted. The GFF file can be sorted using `sort -s -k 1,1 -k 7,7` for strand
    specific, or `sort -s -k 1,1` if not strand specific.

Changes
*******

.. versionchanged:: 0.3.1
    added *cov* command

.. versionchanged:: 0.3.0
    added *--split* option to *sequence* command

.. versionchanged:: 0.2.6
    added *split* command, *--indent* option to *mongodb*

.. versionchanged:: 0.2.3
    added *--gene-id* option to *gtf* command

.. versionadded:: 0.2.2
    added *gtf* command

.. versionadded:: 0.2.1
    *dbm* and *mongodb* commands

.. versionadded:: 0.1.15

"""

from __future__ import division
import sys
import argparse
import logging
import functools
import json

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
        '-s',
        '--split',
        action='store_true',
        help='''Split the sequence header of the reference at the first
        space, to emulate BLAST behaviour''',
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

    seqs = dict(
        (
            seq_id.split(' ')[0] if options.split else seq_id,
            seq
        )
        for seq_id, seq in fasta.load_fasta(options.reference)
    )

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
    parser.add_argument(
        '-i',
        '--indent',
        action='store',
        type=int,
        default=None,
        help='If used, the json will be written in a human readble form'
    )

    parser.set_defaults(func=mongodb_command)


def mongodb_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

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
            annotation.to_mongodb(
                lineage_func=lineage_func,
                indent=options.indent
            )
        )
        options.output_file.write('\n')


def gtf_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    for annotation in gff.parse_gff(options.input_file):
        options.output_file.write(
            annotation.to_gtf(
                gene_id_attr=options.gene_id
            )
        )


def set_gtf_parser(parser):
    parser.add_argument(
        '-g',
        '--gene-id',
        type=str,
        default='gene_id',
        help='GFF attribute to use for the GTF gene_id attribute'
    )
    parser.set_defaults(func=gtf_command)


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


def set_split_parser(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input GFF file, defaults to stdin'
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
        help='Number of chunks into which split the GFF file'
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
        "Splitting GFF into %d chunks with prefix %s",
        options.number,
        options.prefix
    )

    name_mask = "%s-{0:05}.gff" % options.prefix
    if options.gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    gff.split_gff_file(
        options.input_file,
        name_mask,
        num_files=options.number
    )


def set_coverage_parser(parser):
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
        default='-',
        help='Output file (tab separated), defaults to stdout'
    )
    parser.add_argument(
        '-f',
        '--reference',
        type=argparse.FileType('r'),
        required=True,
        help='Reference FASTA file for the GFF'
    )
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        '-j',
        '--json',
        action='store_true',
        default=False,
        help='The output will be a JSON dictionary'
    ),
    group.add_argument(
        '-s',
        '--strand-specific',
        action='store_true',
        default=False,
        help='If the coverage must be calculated on each strand'
    )
    parser.add_argument(
        '-r',
        '--rename',
        action='store_true',
        default=False,
        help='Emulate BLAST output (use only the header part before the' +
        ' first space)'
    )
    parser.set_defaults(func=coverage_command)


def coverage_command(options):
    sequences = dict(
        fasta.load_fasta_rename(options.reference) if options.rename else fasta.load_fasta(options.reference)
    )
    iterator = gff.annotation_coverage_sorted(
        gff.parse_gff(options.input_file),
        sequences,
        strand=options.strand_specific
    )

    contig_coverage = {}

    for seq_id, strand, coverage in iterator:
        if options.json:
            contig_coverage[seq_id] = coverage
        else:
            options.output_file.write(
                "{}\t{}\t{}\n".format(
                    seq_id,
                    "NA" if strand is None else strand,
                    coverage
                )
            )

    if options.json:
        json.dump(contig_coverage, options.output_file, indent=4)


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
    utils.add_basic_options(parser_s, manual=__doc__)

    parser_d = subparsers.add_parser(
        'dbm',
        help='Creates a dbm database'
    )

    set_dbm_parser(parser_d)
    utils.add_basic_options(parser_d, manual=__doc__)

    parser_m = subparsers.add_parser(
        'mongodb',
        help='Extract annotations from a GFF file and makes output for MongoDB'
    )

    set_mongodb_parser(parser_m)
    set_common_options(parser_m)
    utils.add_basic_options(parser_m, manual=__doc__)

    parser_gtf = subparsers.add_parser(
        'gtf',
        help='Extract annotations from a GFF file to a GTF file'
    )

    set_gtf_parser(parser_gtf)
    set_common_options(parser_gtf)
    utils.add_basic_options(parser_gtf, manual=__doc__)

    parser_split = subparsers.add_parser(
        'split',
        help='Split annotations from a GFF file to a several files'
    )

    set_split_parser(parser_split)
    utils.add_basic_options(parser_split, manual=__doc__)

    parser_coverage = subparsers.add_parser(
        'cov',
        help='Report on how much a sequence length is covered by annotations'
    )

    set_coverage_parser(parser_coverage)
    utils.add_basic_options(parser_coverage, manual=__doc__)

    utils.add_basic_options(parser, manual=__doc__)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
