#!/usr/bin/env python
"""
This script downloads sequence data for each gene of interest (some are black
listed because not relevant to the rumen microbiome) and all the specified taxa.

The way a taxon is specified is through a few different rules:

* specific **taxon ids** in uniprot
* a specific **taxon rank** (e.g.: genus, phylum, etc.)
* optional lineage filter: the lineage filter make sure that the name specified
  is included in the lineage attribute in the taxonomy.

As an example, if the rank chosen is genus, and the lineage option is
set to archaea, only the taxa whose rank is genus and that belong to the archaea
subtree will be downloaded:

`download_profiles -m email -r genus -l archaea mg_data/kegg.pickle
mg_data/taxonomy.pickle`

This allows to customise the level of specificity that we want in profiling and
make the process of downloading faster.
For metagenomic data, a good start is mixing different taxon ranks, using the
order or genus for the genes and then specifying a lineage of interest.

Because each profile is indipendent from each other, it's useful to start the
download with a certain rank and then run the profiling. During the profiling
a new download can be started and so on.
"""

import os
import glob
import argparse
import urllib2
import logging
import itertools
import numpy
import pickle
from cStringIO import StringIO
import mgkit
from mgkit.io.fasta import load_fasta
from mgkit import logger
from mgkit import kegg
from mgkit.net import uniprot
from mgkit import taxon
from ..filter.taxon import filter_taxonomy_by_lineage, filter_taxonomy_by_rank
import mgkit.simple_cache
from . import utils


is_ancestor = mgkit.simple_cache.memoize(taxon.is_ancestor)

LOG = logging.getLogger(__name__)


class NotEnoughSeqs(Exception):
    "Raised if only one sequence was downloaded"
    pass


def set_parser():
    "argument parser configuration"
    parser = argparse.ArgumentParser(
        description='Download KO sequences from Uniprot',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        default='profile_files',
        action='store',
        help='directory in which to store the downloaded files'
    )
    parser.add_argument(
        '-k',
        '--kegg-data',
        default='data/kegg.pickle',
        action='store',
        help='pickle file containing Kegg data'
    )
    parser.add_argument(
        '-m',
        '--email',
        required=True,
        action='store',
        help='email address to use for Uniprot communications'
    )
    parser.add_argument(
        '-t',
        '--taxon-data',
        default='data/taxonomy.pickle',
        action='store',
        help='pickle file containing taxonomy data'
    )
    parser.add_argument(
        '-r',
        '--taxon-rank',
        default='genus',
        action='store',
        help='taxon rank to download'
    )
    parser.add_argument(
        '-l',
        '--lineage',
        default='archaea',
        action='store',
        help='lineage for filtering (e.g. archaea)'
    )
    parser.add_argument(
        '-i',
        '--taxon-id',
        action='store',
        nargs='+',
        type=int,
        help='id(s) of taxa to download. If specified take precedence over ' +
        'lineage-rank'
    )
    parser.add_argument(
        '-R',
        '--only-reviewed',
        default=False,
        action='store_true',
        help="Only download reviewed sequences"
    )
    parser.add_argument(
        '-a',
        '--all-path',
        default=False,
        action='store_true',
        help='Download all KO from Kegg - exclude blacklist'
    )
    utils.add_basic_options(parser)

    return parser


def load_data(kegg_data_name, taxon_data, length_data_name):
    "Loads data for script"
    kegg_data = kegg.KeggData(kegg_data_name)
    taxonomy = taxon.UniprotTaxonomy(taxon_data)

    try:
        LOG.info("Using profile's length data (%s)", length_data_name)
        length_data = pickle.load(open(length_data_name, 'r'))
    except IOError:
        LOG.info("Not found (%s), restarting", length_data_name)
        length_data = {}

    return kegg_data, taxonomy, length_data


def choose_taxa(taxonomy, options):
    "Returns the list of ids to look for in Uniprot"
    if options.taxon_id:
        taxon_ids = set(options.taxon_id)
        LOG.debug(taxon_ids)
    else:
        taxon_ids = filter_taxonomy_by_rank(
            taxonomy, options.taxon_rank.lower()
        )
        taxon_ids = filter_taxonomy_by_lineage(
            taxon_ids, options.lineage.lower()
        )
        taxon_ids = set(x.taxon_id for x in taxon_ids)
    LOG.info("The number of taxa to download is %d", len(taxon_ids))
    return taxon_ids


def choose_ko_ids(kegg_data, options):
    "Returns the list of mapping ID->Name according to the options passed"
    if options.all_path:
        LOG.info("Dowloading all KOs, ignoring blacklist")
        ko_ids = kegg_data.get_ko_names()
    else:
        LOG.info("Dowloading only KOs in pathways not blacklisted")
        filt_ko = set(
            itertools.chain(
                *(kegg_data[path_id].genes.keys()
                  for path_id in kegg_data if not path_id in kegg.BLACK_LIST)
            )
        )
        ko_ids = dict(
            (ko_id, name)
            for ko_id, name in kegg_data.get_ko_names().items()
            if ko_id in filt_ko
        )
    return ko_ids


def map_ko_to_uniprot(ko_id, taxon_ids, reviewed, contact):
    "Returns the taxon IDs found in Uniprot for a specific id"
    query_base = "database:(type:ko {0}) AND ({1}){2}"

    query = query_base.format(
        ko_id,
        ' OR '.join(
            "taxonomy:{0}".format(taxon_id) for taxon_id in taxon_ids
        ),
        'AND reviewed:yes' if reviewed else ''
    )

    LOG.info("Downloading taxa list for KO: %s", ko_id)
    taxon_ids = uniprot.ko_to_mapping(
        ko_id, query,
        'organism-id',
        contact=contact
    )

    LOG.info("Found %d taxa for KO: %s", len(taxon_ids), ko_id)

    return taxon_ids


def filter_found_taxa(taxon_ids_found, taxon_ids, taxonomy):
    """
    Filter the taxa found in Uniprot, making sure that they at a lower level
    of those requested
    """
    taxon_ids_found = set(int(taxon_id) for taxon_id in taxon_ids_found)

    taxon_ids_download = set()

    for taxon_id_found in taxon_ids_found:
        if taxon_id_found not in taxonomy:
            LOG.warning(
                "Taxon id (%d) not found in local taxonomy, skipping",
                taxon_id_found
            )
            continue
        for taxon_id in taxon_ids:
            if is_ancestor(taxonomy, taxon_id_found, taxon_id):
                taxon_ids_download.add(taxon_id)
                break

    LOG.info('Downloading a total of %d taxa', len(taxon_ids_download))

    return taxon_ids_download


def download_ko_sequences(ko_id, taxon_ids, reviewed, contact):
    "Downloads the sequences associated to all taxon IDs provided"
    ko_seqs = {}

    for taxon_id in taxon_ids:
        key = (ko_id, taxon_id, reviewed)
        seqs = uniprot.get_sequences_by_ko(
            ko_id,
            taxon_id,
            contact,
            reviewed
        )
        if seqs.count('>') <= 1:
            raise NotEnoughSeqs(
                'Not enough Sequences ({0})'.format(
                    seqs.count('>')
                )
            )
        LOG.debug(
            "Downloaded %d sequences for taxon id: %d",
            seqs.count('>'),
            taxon_id
        )
        ko_seqs[key] = seqs

    return ko_seqs


def write_ko_sequences(seqs, taxonomy, output_dir):
    "Writes fasta sequences to disc"
    for ko_id, taxon_id, reviewed in seqs:
        file_name = '{0}_{1}_{2}{3}.fa'.format(
            ko_id,
            taxon_id,
            taxonomy[taxon_id].s_name.replace(' ', '#'),
            '' if reviewed else '-nr'
        )
        file_handle = open(os.path.join(output_dir, file_name), 'w')
        file_handle.write(seqs[(ko_id, taxon_id, reviewed)])
        file_handle.close()


def add_profiles_to_length(seqs, length_data):
    "Adds the average profile length to the dictionary"
    for key in seqs:
        avg_len = numpy.fromiter(
            (len(seq[1]) for seq in load_fasta(StringIO(seqs[key]))),
            dtype=int
        ).mean()
        length_data[key] = avg_len


def main():
    "Main function"
    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    ret_value = 0

    contact = options.email
    reviewed = options.only_reviewed

    output_dir = options.output_dir
    try:
        LOG.info("Writing fasta files to %s", output_dir)
        os.mkdir(output_dir)
    except OSError:
        LOG.info("Output directory %s already present", output_dir)

    length_data_name = options.output_dir.rstrip('/') + '-length.pickle'

    kegg_data, taxonomy, length_data = load_data(
        options.kegg_data,
        options.taxon_data,
        length_data_name
    )

    taxon_ids = choose_taxa(taxonomy, options)
    ko_ids = choose_ko_ids(kegg_data, options)

    for idx, ko_id in enumerate(sorted(ko_ids)):
        if ko_id in (key[0] for key in length_data):
            LOG.info(
                "(%05d/%05d) Skipped %s, already downloaded",
                idx + 1,
                len(ko_ids),
                ko_id
            )
            continue
        LOG.info("(%05d/%05d) Downloading %s", idx + 1, len(ko_ids), ko_id)
        try:
            taxon_ids_found = map_ko_to_uniprot(ko_id, taxon_ids, reviewed,
                                                contact)
            taxon_ids_found = filter_found_taxa(
                taxon_ids_found,
                taxon_ids,
                taxonomy
            )
            seqs = download_ko_sequences(
                ko_id,
                taxon_ids_found,
                reviewed,
                contact
            )
        except urllib2.HTTPError:
            LOG.error("Couldn't download sequences for %s", ko_id)
            ret_value = 1
            continue
        except NotEnoughSeqs, err:
            LOG.info("%s", err)
            continue

        try:
            write_ko_sequences(seqs, taxonomy, output_dir)
            add_profiles_to_length(seqs, length_data)
        except IOError, error:
            LOG.error("Couldn't write sequences for %s", ko_id)
            LOG.error(str(error))
            file_list = glob.glob(os.path.join(output_dir,
                                  '{0}*'.format(ko_id)))
            for fname in file_list:
                LOG.error("-> %s", fname)
                os.remove(fname)
            ret_value = 2
            continue
        else:
            LOG.debug("Profiles for %s saved", ko_id)
            pickle.dump(length_data, open(length_data_name, 'w'))

    return ret_value

if __name__ == '__main__':
    main()
