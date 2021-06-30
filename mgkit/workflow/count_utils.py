"""
Count Table Utilities

Map Count Table to Genes
************************

The `map` command can map information from map files to create count tables
from featureCounts where `uid` was used as attribute for the counts.

A taxonomy map can be passed if the taxonomy needs to be included in the
index of the output table. The format used for the table is `Parquet`, which
retains the Index/MultiIndex when read back with Pandas

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


def stream_feature_counts(count_file, sample_func=None):
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
        length = int(length)
        counts = pandas.Series([int(count) for count in counts], index=sample_ids, dtype=pandas.SparseDtype('int', 0))
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
    
    iterator = stream_feature_counts(count_file)

    count_data = {}
    pbar = tqdm(desc="Mappings", total=len(map_data))

    for key, length, count in iterator:
        if key not in map_data:
            continue
        for map_id in map_data[key]:
            if taxa_map is not None:
                map_id = (map_id, taxa_data.get(key, None))
            try:
                count_data[map_id] = count_data[map_id] + count
            except KeyError:
                count_data[map_id] = count
        
        # needs to keep it out of the loop
        pbar.update(1)
        
        if pbar.n == len(map_data):
            pbar.close()
            LOG.info("Used all the data from the Map File")
            break
    
    dataframe = pandas.DataFrame.from_dict(count_data, orient='index')
    dataframe.to_parquet(output_file, engine='pyarrow')
