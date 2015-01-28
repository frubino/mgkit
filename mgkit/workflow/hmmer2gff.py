"""
Script to convert HMMER results files (domain table) to a GFF file, the name of
the profiles are expected to be now in the form
*GENEID_TAXONID_TAXON-NAME(-nr)*, but the old one *GENEID_TAXON-NAME(-nr)* is
accepted for compatibility with old experiments.

It is tested only Kegg Orthologs, but the script only uses them if a pickle
file for Kegg is passed to add descriptions. It is general enough to use any
kind of gene identifier.

.. note::

    for GENEID, old documentation points to KOID, it is the same

Changes
*******

.. versionchanged:: 0.1.15
    adapted to new GFF module and specs

"""

import sys
import logging
import argparse
from mgkit import logger
from mgkit.io import gff
from mgkit.io import fasta
from mgkit import kegg
from mgkit.taxon import MISPELLED_TAXA
from . import utils

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Setup command line options
    """
    parser = argparse.ArgumentParser(
        description='Convert HMMER data to GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    utils.add_basic_options(parser)

    group = parser.add_argument_group('File options')
    group.add_argument(
        'aa_file',
        type=argparse.FileType('r'),
        help="Fasta file containing contigs translated to aa (used by HMMER)"
    )
    group.add_argument(
        'hmmer_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-'
    )
    group.add_argument(
        '-o',
        '--output-file',
        nargs='?',
        type=argparse.FileType('w'), default=sys.stdout
    )

    group = parser.add_argument_group('Misc')
    group.add_argument('-q', '--quiet', action='store_const',
                       const=logging.WARNING, default=logging.DEBUG,
                       help='only show warnings or errors')

    group.add_argument('-t', '--discard', action='store', type=float,
                       default=0.05,
                       help='Evalue over which an hit will be discarded')
    group.add_argument('-k', '--ko-descriptions',
                       help='pickle file containing KO descriptions')

    return parser


def get_seq_data(f_handle):
    """
    Load reference sequences.
    """
    # LOG.info('Loading contigs data from file %s', f_handle.name)

    seq_data = dict(
        (name, seq) for name, seq in fasta.load_fasta(f_handle)
    )

    return seq_data


def get_aa_data(f_handle):
    """
    Load aminoacid seuqnces used by HMMER.
    """
    # LOG.info('Loading aa data from file %s', f_handle.name)

    aa_seqs = dict((name, seq) for name, seq in fasta.load_fasta(f_handle))

    return aa_seqs


def parse_domain_table_contigs(f_handle, aa_seqs, f_out, discard,
                               ko_names=None):
    """
    Parse the HMMER result file
    """
    LOG.info('Parsing HMMER data from file %s', f_handle.name)
    LOG.info('Writing GFF data to file %s', f_out.name)

    count_dsc = 0
    count_tot = 0
    count_mis = 0
    count_skp = 0

    ko_counts = {}

    for idx, line in enumerate(f_handle):

        if line.startswith('#'):
            continue
        if idx % 10000 == 0:
            LOG.info("Line number: %d", idx)

        count_tot += 1

        try:
            annotation = gff.from_hmmer(
                line,
                aa_seqs,
                ko_counts=ko_counts
            )
        except ZeroDivisionError:
            LOG.error(
                "Skipping line %d because of an error in the calculations",
                idx + 1
            )
            count_skp += 1
            continue

        if annotation.score > discard:
            count_dsc += 1
            continue
        try:
            annotation.attr['description'] = ko_names[
                annotation.gene_id
            ]
        except (KeyError, TypeError):
            annotation.attr['description'] = ''

        #correct mispelled taxa: profiles4
        try:
            annotation.attr['taxon'] = MISPELLED_TAXA[
                annotation.attr['taxon']
            ]
            LOG.debug("Fixed mispelled taxon %s", annotation.attr['taxon'])
            count_mis += 1
        except KeyError:
            #not mispelled
            pass

        annotation.to_file(f_out)

    LOG.info(
        "Read %d lines, discared %d, skipped %d, mispelled %d",
        count_tot,
        count_dsc,
        count_skp,
        count_mis
    )


def main():
    """
    Main loop
    """
    options = set_parser().parse_args()
    logger.config_log(options.quiet)
    log = logging.getLogger(__name__)
    aa_data = get_aa_data(options.aa_file)

    if options.ko_descriptions:
        log.info("Loading KO descriptions from file %s",
                 options.ko_descriptions)
        kegg_data = kegg.KeggData(options.ko_descriptions)
        ko_names = kegg_data.get_ko_names()
    else:
        ko_names = None

    parse_domain_table_contigs(
        options.hmmer_file,
        aa_data,
        options.output_file,
        options.discard,
        ko_names
    )

if __name__ == '__main__':
    main()
