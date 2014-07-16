"""
Add more information to GFF annotations: gene mappings, coverage, taxonomy,
etc..

Uniprot Command
***************

If the *gene_id* of an annotation is a Uniprot ID, the script queries Uniprot
for the requested information. At the moment the information that can be added
is the taxon_id, taxon_name, lineage and mapping to EC, KO, eggNOG IDs.

It's also possible to add mappings to other databases using the *-m* option with
the correct identifier for the mapping, which can be found at `this page
<http://www.uniprot.org/faq/28>`_; for example if it's we want to add the
mappings of uniprot IDs to *BioCyc*, in the *abbreviation* column of the
mappings we find that it's identifier is *REACTOME_ID*, so we pass *-m REACTOME*
to the script (leaving *_ID* out). Mapped IDs are separated by commas.

The taxonomy IDs are not overwritten if they are found in the annotations, the
*-f* is provided to force the overwriting of those values.

See also :ref:`gff-specs` for more informations about the GFF specifications
used.

.. note::

    As the script needs to query Uniprot a lot, it is reccommended to split
    the GFF in several files, so an error in the connection doesn't waste time.

Taxonomy Command
****************

To refine the taxonomic assignments of predicted genes annotations, the
annotation sequences may be searched against a database like the NCBI *nt*.

This commands takes as input a GFF file, one or more blast output files and a
file with all mappings from GIDs to taxonomy IDs. More information on how to
get the file can be read in the documentation of the function
:func:`mgkit.io.blast.parse_gi_taxa_table`.

The fasta sequences used with BLAST must have as name the uid of the annotations
they refer to, and one way to obtain these sequences is to use the function
:func:`mgkit.io.gff.extract_nuc_seqs` and save them to a fasta file.

The command accept a minimum bitscore to accept an hit and a the best hit is
selected from all those found for a sequence, more details can be found in the
documentation for :func:`mgkit.io.blast.parse_fragment_blast`.

Coverage Command
****************

Adds coverage information from BAM alignment files to a GFF file, using the
function :func:`mgkit.align.add_coverage_info`, the user needs to supply for
each sample a BAM file, using the `-a` option, whose parameter is in the form
`sample,samplealg.bam`. More samples can be supplied adding more `-a` arguments.

.. hint::

    As an example, to add coverage for `sample1`, `sample2` the command line
    is::

        add-gff-info coverage -a sample1,sample1.bam -a sample2,sample2.bam \\
        inputgff outputgff

A total coverage for the annotation is also calculated and stored in the
`cov` attribute, while each sample coverage is stored into `sample_cov` as per
:ref:`gff-specs`.

Changes
*******

.. versionadded:: 0.1.12

.. versionchanged:: 0.1.13
    added *taxonomy* command

.. versionchanged:: 0.1.13
    added *--force-taxon-id* option to the *uniprot*
    command

.. versionchanged:: 0.1.13
    added *coverage* command

"""
from __future__ import division
import sys
import argparse
import logging
import itertools
import pysam
from .. import align
from .. import logger
from . import utils
from ..io import gff, blast
from ..net import uniprot

LOG = logging.getLogger(__name__)


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
        help='Output GFF file, defaults to stdout'
    )


def set_uniprot_parser(parser):
    parser.add_argument(
        '-c',
        '--email',
        action='store',
        type=str,
        help='Contact email',
        default=None
    )
    parser.add_argument(
        '--buffer',
        action='store',
        type=int,
        help='Number of annotations to keep in memory',
        default=50
    )
    group = parser.add_argument_group('Require Internet connection')
    group.add_argument(
        '-f',
        '--force-taxon-id',
        action='store_true',
        default=False,
        help='Overwrite taxon_id if already present'
    )
    group.add_argument(
        '-t',
        '--taxon-id',
        action='store_true',
        default=False,
        help="""
             Add taxonomic ids to annotations, if taxon_id is found, it won't
             be Overwritten.
             """
    )
    group.add_argument(
        '-l',
        '--lineage',
        action='store_true',
        default=False,
        help='Add taxonomic lineage to annotations'
    )
    group.add_argument(
        '-e',
        '--eggnog',
        action='store_true',
        default=False,
        help='Add eggNOG mappings to annotations'
    )
    group.add_argument(
        '-ec',
        action='store_true',
        default=False,
        help='Add EC mappings to annotations'
    )
    group.add_argument(
        '-ko',
        action='store_true',
        default=False,
        help='Add KO mappings to annotations'
    )
    group.add_argument(
        '-m',
        '--mapping',
        action='append',
        type=str,
        help='Add any DB mappings to annotations'
    )

    parser.set_defaults(func=uniprot_command)


def add_uniprot_info(annotations, options):
    columns = []
    if options.taxon_id:
        columns.append('organism')
        columns.append('organism-id')
    if options.eggnog:
        columns.append('database(EGGNOG)')
    if options.ko:
        columns.append('database(KO)')
    if options.ec:
        columns.append('ec')
    if options.lineage:
        columns.append('lineage()')
    if options.mapping is not None:
        for db in options.mapping:
            columns.append('database({0})'.format(db))

    if not columns:
        return

    LOG.info("Retrieving gene information from Uniprot")

    data = uniprot.get_gene_info(
        [x.gene_id for x in annotations],
        columns=columns,
        contact=options.email
    )

    for annotation in annotations:
        try:
            gene_info = data[annotation.gene_id]
        except KeyError:
            #no data was found
            continue

        for column, values in gene_info.iteritems():
            #nothing found
            if not values:
                continue

            if column == 'organism-id':
                if (annotation.taxon_id and options.force_taxon_id) or \
                   (annotation.taxon_id is None):
                    annotation.attr['taxon_id'] = int(values)
                    annotation.attr['taxon_db'] = 'UNIPROT'
                    #test with a try/expect maybe
                    if 'organism' in columns:
                        annotation.attr['taxon_name'] = gene_info['organism']
                    if column.startswith('lineage'):
                        annotation.attr['lineage'] = gene_info['lineage()']
            elif column.startswith('database'):
                annotation.attr[
                    'map_{0}'.format(
                        column[:-1].split('(')[1]
                    )
                ] = ','.join(values)
            elif column == 'ec':
                if isinstance(values, list):
                    annotation.attr['EC'] = ','.join(
                        x.split('-')[0]
                        for x in values
                    )
                else:
                    annotation.attr['EC'] = values


def uniprot_command(options):

    if options.buffer < 1:
        options.buffer = 1

    ann_buffer = []

    for annotation in gff.parse_gff(options.input_file, gff_type=gff.from_gff):

        ann_buffer.append(annotation)

        if len(ann_buffer) == options.buffer:

            add_uniprot_info(ann_buffer, options)

            for annotation in ann_buffer:
                annotation.to_file(options.output_file)

            ann_buffer = []
    else:
        add_uniprot_info(ann_buffer, options)

        for annotation in ann_buffer:
            annotation.to_file(options.output_file)


def set_blast_taxonomy_parser(parser):
    parser.add_argument(
        '-t',
        '--gi-taxa-table',
        action='store',
        default=None,
        required=True,
        help="GIDs taxonomy table (e.g. gi_taxid_nucl.dmp.gz)"
    )
    parser.add_argument(
        '-b',
        '--blast-output',
        action='store',
        nargs='+',
        required=True,
        help="BLAST output file(s)"
    )
    parser.add_argument(
        '-s',
        '--bitscore',
        action='store',
        default=40,
        type=float,
        help="Minimum bitscore allowed"
    )
    parser.add_argument(
        '-d',
        '--taxon-db',
        action='store',
        default='NCBI-NT',
        help="NCBI database used"
    )
    parser.set_defaults(func=taxonomy_command)


def taxonomy_command(options):

    uid_gid_map = dict(
        itertools.chain(
            *(blast.parse_fragment_blast(x, bitscore=options.bitscore)
              for x in options.blast_output)
        )
    )

    gids = set(x[0] for x in uid_gid_map.itervalues())

    gid_taxon_map = dict(
        (gid, taxon_id)
        for gid, taxon_id in blast.parse_gi_taxa_table(
            options.gi_taxa_table, gids=gids
        )
    )

    count = 0
    tot_count = 0

    for annotation in gff.parse_gff(options.input_file):
        tot_count += 1
        try:
            gid, bitscore, identity = uid_gid_map[annotation.uid]
            annotation.taxon_id = gid_taxon_map[gid]
            annotation.taxon_db = options.taxon_db
        except KeyError:
            continue
        finally:
            annotation.to_file(options.output_file)

        count += 1

    LOG.info(
        "Added taxonomy information for %.2f%% annotations (%d/%d)",
        count / tot_count * 100,
        count,
        tot_count
    )


def split_sample_alg(argument):
    "Split sample/alignment option"
    try:
        sample, bam_file_name = argument.split(',', 1)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Can't get get both sample and bam file from '%s'" % argument
        )

    return (sample, bam_file_name)


def set_coverage_parser(parser):
    parser.add_argument(
        '-a',
        '--sample-alignment',
        action='append',
        required=True,
        type=split_sample_alg,
        help='sample name and correspondent alignment file separated by comma'
    )
    parser.set_defaults(func=coverage_command)


def coverage_command(options):
    samples = []
    bam_files = []

    for sample, bam_file_name in options.sample_alignment:
        samples.append(sample)
        bam_files.append(pysam.Samfile(bam_file_name, 'rb'))

    if len(samples) != len(set(samples)):
        LOG.critical("There are duplicate sample names")
        return 1

    annotations = list(gff.parse_gff(options.input_file))

    align.add_coverage_info(annotations, bam_files, samples)

    gff.write_gff(annotations, options.output_file)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Adds informations to a GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_u = subparsers.add_parser(
        'uniprot',
        help='Adds information from GFF whose gene_id is from Uniprot'
    )

    set_uniprot_parser(parser_u)
    set_common_options(parser_u)

    parser_t = subparsers.add_parser(
        'taxonomy',
        help='''Adds taxonomic information from annotation sequences blasted
                against a NCBI db'''
    )

    set_blast_taxonomy_parser(parser_t)
    set_common_options(parser_t)

    parser_c = subparsers.add_parser(
        'coverage',
        help='Adds coverage information from BAM Alignment files'
    )

    set_coverage_parser(parser_c)
    set_common_options(parser_c)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)
