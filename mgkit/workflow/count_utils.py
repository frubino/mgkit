"""
Count Table Utilities

Map Count Table to Genes
************************

The `map` command can map information from map files to create count tables
from featureCounts where `uid` was used as attribute for the counts.

A taxonomy map can be passed if the taxonomy needs to be included in the
index of the output table. The format used for the table is `Parquet`, which
retains the Index/MultiIndex when read back with Pandas

Concatenate Parquet Files
*************************

Allows to concatenate several pandas dataframe with same indices. It's used
when the mapping file produce too big files and won't fit in memory.

So a solution is to split the map files, making multiple parquet files and
after that, concatenate them with this script.

Convert Parquet into CSV
************************

The command `to_csv` outputs a CSV file from a Parquet table.

Changes
*******

.. versionadded:: 0.5.7

"""

import logging
import click
from . import utils
from .. import logger
from tqdm import tqdm
import pandas
from mgkit.io import gff, open_file
from mgkit.utils.dictionary import text_to_dict


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"


def stream_feature_counts(count_file, sample_func=None, gene_ids=None):
    LOG.info("Reading featureCounts file %s", count_file)
    # comment at the start
    _ = count_file.readline()
    # header
    _, _, _, _, _, _, *sample_ids = count_file.readline().strip().split('\t')
    if sample_func is not None:
        sample_ids = [sample_func(sample_id) for sample_id in sample_ids]
    
    for line in count_file:
        if line.startswith('#'):
            continue
        gene_id, seq_id, start, end, strand, length, *counts = line.strip().split('\t')
        # skip rows not in the map
        if (gene_ids is not None) and (gene_id not in gene_ids):
            continue
        length = int(length)
        counts = pandas.Series((int(count) for count in counts), index=sample_ids, dtype=pandas.SparseDtype('int', 0))
        yield gene_id, length, counts


@main.command('map', help="""Map counts with information a dictionary file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-m', '--map-file', type=click.File('r'), required=True,
                help="Map file to use")
@click.option('-t', '--taxa-map', type=click.File('r'), default=None,
                help="Taxa map file")
@click.option('-s', '--separator', default="\t", show_default=True, type=click.STRING,
                help="Field separator for map file Key/Value")
@click.option('-sv', '--split-value', default=False, show_default=True, is_flag=True,
              help="Values are string to be split")
@click.argument('count_file', type=click.File('r', lazy=False), default='-')
@click.argument('output_file', type=click.File('wb'), default='-')
def map_counts(verbose, map_file, taxa_map, separator, split_value, count_file, output_file):
    """
    .. versionadded: 0.5.7

    Map counts from a key specified in the featureCounts file (first column) to another
    """
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if taxa_map is not None:
        taxa_data = {}
        LOG.info("Reading Taxa Map file %s", taxa_map)
        for line in taxa_map:
            key, taxon = line.strip().split(separator)
            taxa_data[key] = taxon
        LOG.info("Finished reading Taxa Map file")

    map_data = {}
    LOG.info("Reading Map file %s", map_file)
    for line in map_file:
        key, *values = line.strip().split(separator)
        if split_value:
            values = list(values[0])
        map_data[key] = set(values)
    LOG.info("Finished reading Map file")
    
    iterator = stream_feature_counts(count_file, gene_ids=map_data)

    count_data = {}
    pbar = tqdm(desc="Mappings", total=len(map_data))

    for key, length, count in iterator:
        # keys was found, so needs to update
        pbar.update(1)

        for map_id in map_data[key]:
            if taxa_map is not None:
                map_id = (map_id, taxa_data.get(key, None))
            try:
                count_data[map_id] = count_data[map_id] + count
            except KeyError:
                count_data[map_id] = count
        
        # If all keys were found, break the loop
        if pbar.n == len(map_data):
            pbar.close()
            LOG.info("Used all the data from the Map File")
            break

    LOG.info("Building DataFrame")
    count_data = pandas.DataFrame.from_dict(count_data, orient='index').sort_index()
    count_data.to_parquet(output_file, engine='pyarrow')


@main.command('cat', help="""Combine multiple count tables files""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-o', '--output', required=True, help="Output file",
                type=click.Path(file_okay=True, writable=True))
@click.argument('count_files', type=click.Path(file_okay=True, readable=True), nargs = -1)
def cat_count_files(verbose, output, count_files):
    """
    .. versionadded: 0.5.7

    Concatenate multiple count tables - use when there are memory constrains in bigger count
    tables
    """

    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    dataframes = []
    for count_file in count_files:
        LOG.info("Reading file %s", count_file)
        dataframes.append(pandas.read_parquet(count_file))

    n_levels = dataframes[0].index.nlevels
    keys = list(range(len(dataframes)))

    dataframe = pandas.concat(
        dataframes,
        keys=keys,
    ).groupby(level=list(range(1, n_levels + 1))).sum()

    dataframe.to_parquet(output)

@main.command('to_csv', help="Convert Parquet tables into CSV")
@click.option('-v', '--verbose', is_flag=True)
@click.argument('parquet_file', type=click.Path(file_okay=True, readable=True))
@click.argument('csv_file', type=click.File('w'), default='-')
def to_csv_command(verbose, parquet_file, csv_file):
    """
    .. versionadded: 0.5.7

    Converts a Parquet table to CSV
    """
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    dataframe = pandas.read_parquet(parquet_file)
    dataframe.to_csv(csv_file)
