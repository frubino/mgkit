"""
Annotate GFF with more information/mappings
"""
from __future__ import division
import sys
import argparse
import logging
import itertools
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

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)
