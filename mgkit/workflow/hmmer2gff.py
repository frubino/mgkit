"""
Script to convert HMMER results files (domain table) to a GFF file, the name of
the profiles are expected to be now in the form
*GENEID_TAXONID_TAXON-NAME(-nr)* by default, but any other profile name is
accepted.

The profiles tested are those made from Kegg Orthologs, from the
`download_profiles` script. If the `--no-custom-profiles` options is used,
the script can be used with any profile name. The profile name will be used
for `gene_id`, `taxon_id` and `taxon_name` in the GFF file.

It is possible to use seuqnces not translated using mgkit, no information on
the frame is assumed, so this script can be used against a protein DB. For
example Uniprot can be searched for profiles, in which case the **--no-frame**
options must be used.

.. note::

    for GENEID, old documentation points to KOID, it is the same

.. warning::

    The compatibility with old data has been **removed**, meaning that old
    experiments must use the scripts from those versions. It is possible to use
    multiple environments, with `virtualenv` for this purpose. An examples is
    given in :ref:`install-ref`.

Changes
*******

.. versionchanged:: 0.1.15
    adapted to new GFF module and specs

.. versionchanged:: 0.2.1
    added options to customise output and filters and old restrictions

.. versionchanged:: 0.3.1
    added *--no-frame* option for non mgkit-translated proteins, sequence
    headers are handled the same way as HMMER (truncated at the first space)

"""

import sys
import logging
import argparse
from mgkit import logger
from mgkit.io import gff
from mgkit.io import fasta
from . import utils

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Setup command line options
    """
    parser = argparse.ArgumentParser(
        description='Convert HMMER data to GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

    group = parser.add_argument_group('Filters')
    group.add_argument(
        '-t',
        '--discard',
        action='store',
        type=float,
        default=0.05,
        help='Evalue over which an hit will be discarded'
    )
    group.add_argument(
        '-d',
        '--disable-evalue',
        action='store_true',
        default=False,
        help='Disable Evalue filter'
    )

    group = parser.add_argument_group('GFF')
    group.add_argument(
        '-c',
        '--no-custom-profiles',
        action='store_false',
        default=True,
        help='Profiles names are not in the custom format'
    )
    group.add_argument(
        '-db',
        '--database',
        action='store',
        default='CUSTOM',
        help='Database from which the profiles are generated " +" (e.g. PFAM)'
    )
    group.add_argument(
        '-f',
        '--feature-type',
        action='store',
        default='gene',
        help='Type of feature (e.g. gene)'
    )
    group.add_argument(
        '-n',
        '--no-frame',
        action='store_true',
        default=False,
        help='Set if the sequences were not translated with translate_seq'
    )

    group = parser.add_argument_group('Misc')

    utils.add_basic_options(parser, manual=__doc__)

    return parser


def get_aa_data(f_handle):
    """
    Load aminoacid seuqnces used by HMMER.
    """
    # LOG.info('Loading aa data from file %s', f_handle.name)

    aa_seqs = dict(
        (name.split(' ')[0], seq)
        for name, seq in fasta.load_fasta(f_handle)
    )

    return aa_seqs


def parse_domain_table_contigs(options):
    """
    Parse the HMMER result file
    """
    aa_seqs = get_aa_data(options.aa_file)

    LOG.info('Parsing HMMER data from file %s', options.hmmer_file.name)
    LOG.info('Writing GFF data to file %s', options.output_file.name)

    count_dsc = 0
    count_tot = 0
    count_skp = 0

    for idx, line in enumerate(options.hmmer_file):

        if line.startswith('#'):
            continue
        if idx % 10000 == 0:
            LOG.info("Line number: %d", idx)

        count_tot += 1

        try:
            annotation = gff.from_hmmer(
                line,
                aa_seqs,
                feat_type=options.feature_type,
                db=options.database,
                custom_profiles=options.no_custom_profiles,
                noframe=options.no_frame,
            )
        except ZeroDivisionError:
            LOG.error(
                "Skipping line %d because of an error in the calculations",
                idx + 1
            )
            count_skp += 1
            continue

        # if disable_evalue is True, skips filter
        if not options.disable_evalue:
            if annotation.score > options.discard:
                count_dsc += 1
                continue

        annotation.to_file(options.output_file)

    LOG.info(
        "Read %d lines, discarded %d, skipped %d",
        count_tot,
        count_dsc,
        count_skp
    )


def main():
    """
    Main loop
    """
    options = set_parser().parse_args()
    logger.config_log(options.verbose)

    parse_domain_table_contigs(
        options
    )


if __name__ == '__main__':
    main()
