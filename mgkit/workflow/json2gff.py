"""

.. versionchanged:: 0.3.4
    using *click* instead of *argparse*

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
import click
import logging
import json

import mgkit
from . import utils
from mgkit.io import gff

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('mongodb', help="""Convert annotations from a MongoDB instance to
              GFF""")
@click.option('-v', '--verbose', is_flag=True)
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('gff-file', type=click.File('wb'), default='-')
def mongodb_command(verbose, input_file, gff_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(gff_file, 'name', repr(gff_file))
    )

    LOG.info(
        'Reading from file (%s)',
        getattr(input_file, 'name', repr(input_file))
    )

    for line in input_file:
        line = line.decode('ascii')
        annotation = gff.from_mongodb(json.loads(line), lineage=False)
        annotation.to_file(gff_file)
