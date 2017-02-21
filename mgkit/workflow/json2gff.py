"""
.. versionadded:: 0.2.6

This script converts annotations in JSON format that were created using MGKit
back into GFF annotations.

mongodb command
***************

Annotations converted into MongoDB records with *get-gff-info mongodb* can be
converted back into a GFF file using this command. It can be useful to get a
GFF file as output from a query to a MongoDB instance on the command line.

For example:

mongoexport -d db -c test | json2gff mongodb

will convert all the annotations in the database *db*, collection *test* to
the standard out.

"""

from __future__ import division
import argparse
import logging
import json
import sys

import mgkit
from . import utils
from mgkit.io import gff

LOG = logging.getLogger(__name__)


def set_common_options(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Output file, defaults to stdout'
    )


def set_mongodb_parser(parser):
    parser.set_defaults(func=mongodb_command)


def mongodb_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    LOG.info(
        'Reading from file (%s)',
        getattr(options.input_file, 'name', repr(options.input_file))
    )

    for line in options.input_file:
        annotation = gff.from_mongodb(json.loads(line), lineage=False)
        annotation.to_file(options.output_file)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Convert JSON to GFF',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_mongodb = subparsers.add_parser(
        'mongodb',
        help='Convert annotations from a MongoDB instance to GFF'
    )
    set_mongodb_parser(parser_mongodb)
    set_common_options(parser_mongodb)
    utils.add_basic_options(parser_mongodb, manual=__doc__)

    parser_json = subparsers.add_parser(
        'json',
        help='Convert annotations from JSON lines to GFF'
    )
    set_mongodb_parser(parser_json)
    set_common_options(parser_json)
    utils.add_basic_options(parser_json, manual=__doc__)

    utils.add_basic_options(parser, manual=__doc__)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
