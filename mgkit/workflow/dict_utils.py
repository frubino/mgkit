"""
Dictionary Utilities
********************

This scripts include a series of commands to help with dictionary files, where
the first columns is a key and all columns after are values. Most of the command
expect the structure of a dictionary file like this::

key1,value1,value2,..
key2,value3,value4,..

So all columns after the first are supposed to be part of a list. The internal
representation is a dictionary with keys and values as lists (or sets).

Split
*****

The `split` command can be used in those cases where the value is just one column
and the values must be separated. For example, if you use edit-gff view to extract
the `uid` and `map_KO` in a GFF from MGKit, you'll have::

uid1 K00001,K00002
uid2 K00004,K01002

By passing the value separator, they can be split into the form accepted by other
commands. The same can be accomplished by using `tr` in linux to substitute commas
with tabs.

A trickier and more compelling use is for the `FC` attribute. They are in the form::

uid1 AK
uid2 FU

These categories are single letter and they can split by using the `--no-separator`
option.

Reverse Dictionary
******************

Used to reverse a dictionary Key/Values. Values are used as keys in the new
dictionary and the original keys become the values.

Changes
*******

.. versionadded:: 0.5.7

"""
import logging
import click
import random
from tqdm import tqdm
from . import utils
from .. import logger
from mgkit.utils.dictionary import text_to_dict


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"


@main.command('split', help="""Split values in a dictionary file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--separator', default="\t", show_default=True, type=click.STRING,
                help="Field separator for map file Key/Value")
@click.option('-vs', '--value-separator', default=",", show_default=True, type=click.STRING,
                help="Field separator values")
@click.option('--no-separator', default=False, show_default=True, is_flag=True,
              help="Values are string to be split by character")
@click.option('-os', '--output-separator', default="\t", show_default=True, type=click.STRING,
                help="Field separator for Output map file Key/Value")
@click.argument('input_file', type=click.File('r', lazy=False), default='-')
@click.argument('output_file', type=click.File('w'), default='-')
def group_dict(verbose, separator, value_separator, no_separator, output_separator, 
               input_file, output_file):
    
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if no_separator:
        value_func = list
    else:
        value_func = lambda x: x.split(value_separator)
    
    input_iterator = text_to_dict(input_file, sep=separator, value_func=value_func, 
                                  skip_comment='#', skip_empty=True, verbose=True)

    for key, value in tqdm(input_iterator, desc='Reading File'):
        print(key, *value, sep=output_separator, file=output_file)


@main.command('reverse', help="""Reverse Key/Value in a dictionary file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--separator', default="\t", show_default=True, type=click.STRING,
                help="Field separator for map file Key/Value")
@click.option('-os', '--output-separator', default="\t", show_default=True, type=click.STRING,
                help="Field separator for Output map file Key/Value")
@click.option('-r', '--randomise', is_flag=True)
@click.argument('input_file', type=click.File('r', lazy=False), default='-')
@click.argument('output_file', type=click.File('w'), default='-')
def reverse_dict(verbose, separator, output_separator, randomise, input_file, output_file):

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    reverse_map = {}

    for line in tqdm(input_file, desc='Reading File'):
        key, *values = line.strip().split(separator)
        for value in values:
            try:
                reverse_map[value].add(key)
            except KeyError:
                reverse_map[value] = set([key])
    
    key_list = list(reverse_map.keys())
    if randomise:
        LOG.info('Randomising Keys')
        random.shuffle(key_list)

    for key in tqdm(key_list, desc='Writing File'):
        values = reverse_map[key]
        print(key, *values, sep=output_separator, file=output_file)


