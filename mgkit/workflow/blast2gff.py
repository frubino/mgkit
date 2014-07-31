"""
Blast output conversion in GFF requires a BLAST+ tabular format which can be
obtained by using the `--outfmt 6` option with the default columns, as
specified in :func:`mgkit.io.blast.parse_blast_tab`. The script can get data
from the standard in and ouputs GFF lines on the standard output by default.

Uniprot
*******

The Function :func:`mgkit.io.blast.parse_uniprot_blast` is used, which filters
BLAST hits based on bitscore and adds by default a *db* attribute to the
annotation with the value `UNIPROT-SP`, indicating that the SwissProt db is
used and a *dbq* attribute with the value 10. The feature type used in the GFF
is CDS.

.. blockdiag::

    {
        "BLAST+" [color = "#377eb8" , textcolor = 'white', shape = flowchart.input];
        "parse_uniprot_blast" [color = "#e41a1c" , textcolor = 'white', width=200, fontsize=16];
        "GFF" [color = "#4daf4a" , textcolor = 'white'];
        "BLAST+"  -> "parse_uniprot_blast" -> GFF;
    }

Changes
*******

.. versionadded:: 0.1.12

.. versionchanged:: 0.1.13

* added *-n* parameter to *uniprot* command
* added *-k* option to *uniprot* command

"""
import sys
import argparse
import logging
from .. import logger
from . import utils
from ..io import blast


LOG = logging.getLogger(__name__)


def parse_attr_arg(value):
    values = value.split(':')
    if len(values) != 2:
        raise argparse.ArgumentTypeError(
            "Wrong filter format, must be 'key:value' 'key:value'"
        )

    return values[0], values[1]


def set_common_options(parser):
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
        '-k',
        '--attr-value',
        action='append',
        type=parse_attr_arg,
        help='''Additional attribute and value to add to each annotation,
                in the form attr:value''',
        default=None
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='BLAST+ output file in tabular format, defaults to stdin'
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
        '-db',
        '--db-used',
        action='store',
        type=str,
        default='UNIPROT-SP',
        help='Uniprot database used with BLAST'
    )
    parser.add_argument(
        '-n',
        '--no-split',
        action='store_true',
        default=False,
        help='''if used, the script assumes that the sequence header contains
                only the gene id'''
    )

    parser.set_defaults(func=convert_from_uniprot)


def convert_from_uniprot(options):

    if options.no_split is True:
        name_func = lambda x: x
    else:
        name_func = None

    iterator = blast.parse_uniprot_blast(
        options.input_file,
        bitscore=options.bitscore,
        db=options.db_used,
        dbq=options.db_quality,
        name_func=name_func
    )

    for annotation in iterator:
        if options.attr_value is not None:
            for key, value in options.attr_value:
                annotation.set_attr(key, value)

        annotation.to_file(options.output_file)


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
        help='Blast results from a Uniprot database, by default SwissProt'
    )

    set_uniprot_parser(parser_u)
    set_common_options(parser_u)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)
