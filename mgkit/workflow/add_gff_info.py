"""
Annotate GFF with more information/mappings
"""
import sys
import argparse
import logging
from .. import logger
from . import utils
from ..io import gff
from ..net import uniprot

LOG = logging.getLogger(__name__)


def set_common_options(parser):
    parser.add_argument(
        '-f',
        '--buffer',
        action='store',
        type=int,
        help='Number of annotations to keep in memory',
        default=50
    )
    parser.add_argument(
        '-c',
        '--email',
        action='store',
        type=str,
        help='Contact email',
        default=None
    )
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
    group = parser.add_argument_group('Require Internet connection')
    group.add_argument(
        '-t',
        '--taxon-id',
        action='store_true',
        default=False,
        help='Add taxonomic ids to annotations'
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
            #no data were found
            continue

        for column, values in gene_info.iteritems():
            #nothing found
            if not values:
                continue
            if column == 'organism':
                annotation.attr['taxon_name'] = values
            elif column == 'organism-id':
                annotation.attr['taxon_id'] = int(values)
                annotation.attr['taxon_db'] = 'UNIPROT'
            elif column.startswith('database'):
                annotation.attr[
                    column[:-1].split('(')[1]
                ] = ','.join(values)
            elif column == 'ec':
                if isinstance(values, list):
                    annotation.attr['EC'] = ','.join(
                        x.split('-')[0]
                        for x in values
                    )
                else:
                    annotation.attr['EC'] = values.split('-')[0]
            elif column.startswith('lineage'):
                annotation.attr['lineage'] = values


def uniprot_command(options):

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

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    if options.buffer < 1:
        options.buffer = 1

    logger.config_log(options.verbose)
    options.func(options)
