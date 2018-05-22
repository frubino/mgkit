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

.. versionchanged:: 0.3.4
    using *click* instead of *argparse*, renamed *split* command *--json* to
    *--json-out*

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
import logging
import functools
import json
import click

import mgkit
from . import utils
from mgkit.io import gff, fasta
from mgkit.db import dbm
from mgkit import taxon
from mgkit import simple_cache

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('sequence', help='''Extract the nucleotidic sequences of
annotations from [gff-file] to [fasta-file]''')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-r', '--reverse', is_flag=True,
              help='Reverse complement sequences on the - strand')
@click.option('-w', '--no-wrap', is_flag=True,
              help='Write the sequences on one line')
@click.option('-s', '--split', is_flag=True,
              help='''Split the sequence header of the reference at the first space, to emulate BLAST behaviour''')
@click.option('-f', '--reference', type=click.File('rb'), default=None,
              help='Fasta file containing the reference sequences of the GFF file')
@click.argument('gff-file', type=click.File('rb'), default='-')
@click.argument('fasta-file', type=click.File('wb'), default='-')
def sequence_command(verbose, reverse, no_wrap, split, reference, gff_file,
                     fasta_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if reference is None:
        utils.exit_script('A fasta reference file is required', 1)

    wrap = 60

    if no_wrap:
        wrap = None

    seqs = dict(
        (
            seq_id.split(' ')[0] if split else seq_id,
            seq
        )
        for seq_id, seq in fasta.load_fasta(reference)
    )

    ann_iter = gff.parse_gff(gff_file, gff_type=gff.from_gff)

    seq_iter = gff.extract_nuc_seqs(ann_iter, seqs, reverse=reverse)

    for name, seq in seq_iter:
        fasta.write_fasta_sequence(fasta_file, name, seq, wrap=wrap)


@main.command('dbm', help='''Creates a dbm database with annotations from file
[gff-file] into db [output-dir]''')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-d', '--output-dir', default='gff-dbm', show_default=True,
              help='Directory for the database')
@click.argument('gff-file', type=click.File('rb'), default='-')
def dbm_command(verbose, output_dir, gff_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    db = dbm.create_gff_dbm(gff.parse_gff(gff_file), output_dir)
    db.close()


@main.command('mongodb', help='''Extract annotations from a GFF [gff-file] file
and makes output for MongoDB [output-file]''')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--taxonomy', type=click.File('rb'), default=None,
              help='Taxonomy used to populate the lineage')
@click.option('-c', '--no-cache', is_flag=True,
              help='No cache for the lineage function')
@click.option('-i', '--indent', type=click.INT, default=None,
              help='If used, the json will be written in a human readble form')
@click.argument('gff-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def mongodb_command(verbose, taxonomy, no_cache, indent, gff_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    lineage_func = None

    if taxonomy is not None:
        taxonomy = taxon.Taxonomy(taxonomy)
        lineage_func = functools.partial(
            taxon.get_lineage,
            taxonomy
        )
        if no_cache:
            LOG.info('Using cached calls to lineage')
            lineage_func = simple_cache.memoize(lineage_func)

    for annotation in gff.parse_gff(gff_file):
        output_file.write(
            annotation.to_mongodb(
                lineage_func=lineage_func,
                indent=indent
            ).encode('ascii')
        )
        output_file.write('\n'.encode('ascii'))


@main.command('gtf', help='''Extract annotations from a GFF file [gff-file] to
a GTF file [gtf-file]''')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-g', '--gene-id', default='gene_id', show_default=True,
              help='GFF attribute to use for the GTF *gene_id* attribute')
@click.argument('gff-file', type=click.File('rb'), default='-')
@click.argument('gtf-file', type=click.File('wb'), default='-')
def gtf_command(verbose, gene_id, gff_file, gtf_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(gtf_file, 'name', repr(gtf_file))
    )

    for annotation in gff.parse_gff(gff_file):
        gtf_file.write(
            annotation.to_gtf(
                gene_id_attr=gene_id
            ).encode('ascii')
        )


@main.command('split', help="""Split annotations from a GFF file [gff-file] to
several files starting with [prefix]""")
@click.option('-v','--verbose', is_flag=True,)
@click.option('-p', '--prefix', default='split', show_default=True,
              help='Prefix for the file name in output')
@click.option('-n', '--number', type=click.INT, default=10, show_default=True,
              help='Number of chunks into which split the GFF file')
@click.option('-z', '--gzip', is_flag=True, default=False,
              help='gzip output files')
@click.argument('gff-file', type=click.File('rb'), default='-')
def split_command(verbose, prefix, number, gzip, gff_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        "Splitting GFF into %d chunks with prefix %s",
        number,
        prefix
    )

    name_mask = "%s-{0:05}.gff" % prefix
    if gzip:
        name_mask += '.gz'
        LOG.info("Output files will be compressed (gzip)")

    gff.split_gff_file(
        gff_file,
        name_mask,
        num_files=number
    )


@main.command('cov', help="""Report on how much a sequence length is covered
by annotations in [gff-file]""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-f', '--reference', type=click.File('rb'), required=True,
              help='Reference FASTA file for the GFF')
@click.option('-j', '--json-out', is_flag=True, default=False,
              help='The output will be a JSON dictionary')
@click.option('-s', '--strand-specific', is_flag=True, default=False,
              help='If the coverage must be calculated on each strand')
@click.option('-r', '--rename', default=False, is_flag=True,
              help='Emulate BLAST output (use only the header part before the first space)')
@click.argument('gff-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def coverage_command(verbose, reference, json_out, strand_specific, rename,
                     gff_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    sequences = dict(
        fasta.load_fasta_rename(reference) if rename else fasta.load_fasta(reference)
    )
    iterator = gff.annotation_coverage_sorted(
        gff.parse_gff(gff_file),
        sequences,
        strand=strand_specific
    )

    contig_coverage = {}

    for seq_id, strand, coverage in iterator:
        if json_out:
            contig_coverage[seq_id] = coverage
        else:
            output_file.write(
                "{}\t{}\t{}\n".format(
                    seq_id,
                    "NA" if strand is None else strand,
                    coverage
                ).encode('ascii')
            )

    if json_out:
        output_file.write(
            json.dumps(contig_coverage, indent=4).encode('ascii')
        )
