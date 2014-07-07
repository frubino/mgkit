"""
.. versionadded:: 0.1.12

Blast output conversion in GFF
"""
import sys
import argparse
import logging
from .. import logger
from . import utils
from ..io import blast, gff
from ..net import uniprot

LOG = logging.getLogger(__name__)


def set_common_options(parser):
    parser.add_argument(
        '-bf',
        '--buffer',
        action='store',
        type=int,
        help='Number of annotations to keep in memory',
        default=50
    )
    parser.add_argument(
        '-dbq',
        '--db-quality',
        action='store',
        type=int,
        help='Quality of the DB used',
        default=10
    )
    parser.add_argument(
        '-b',
        '--bitscore',
        action='store',
        type=float,
        help='Minimum bitscore to keep the annotation',
        default=0.0
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
        default='-'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout
    )


def set_uniprot_parser(parser):
    parser.add_argument(
        '-db',
        '--db-used',
        action='store',
        type=str,
        default='UNIPROT-SP',
        help='Uniprot database used with BLAST'
    )
    group = parser.add_argument_group('Requires Internet connection')
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
    parser.set_defaults(func=convert_from_uniprot)


def convert_from_uniprot(options):

    iterator = blast.parse_uniprot_blast(
        options.input_file,
        bitscore=options.bitscore,
        db=options.db_used,
        dbq=options.db_quality
    )

    ann_buffer = []

    for annotation in iterator:
        ann_buffer.append(annotation)

        #write to disk
        if len(ann_buffer) == options.buffer:
            add_uniprot_info(ann_buffer, options)
            gff.write_gff(ann_buffer, options.output_file, verbose=False)
            ann_buffer = []
    else:
        add_uniprot_info(ann_buffer, options)
        gff.write_gff(ann_buffer, options.output_file, verbose=False)


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


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Convert BLAST output to a GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()
    parser_u = subparsers.add_parser(
        'uniprot',
        help='Blast results from Uniprot'
    )

    set_uniprot_parser(parser_u)
    set_common_options(parser_u)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    if options.buffer == 0:
        options.buffer = 1

    logger.config_log(options.verbose)
    options.func(options)
