"""
.. versionadded:: 0.2.5

The script contains commands used to access functionality related to
taxonomy, without the need to write ad-hoc code for functionality that
can be part of a workflow. One example is access to the the last common
ancestor function contained in the :mod:`mgkit.taxon`.

lca and lca_line commands
*************************

These commands expose the functionality of
:func:`last_common_ancestor_multiple`, making it accessible via the command
line. They differ in the input file format and the choice of output files.

the *lca* command can be used to define the last common ancestor of contigs
from the annotation in a GFF file. The command uses the *taxon_ids* from all
annotations belonging to a contig/sequence, if they have a **bitscore** higher
or equal to the one passed (50 by default). The default output of the command
is a tab separated file where the first column is the contig/sequence name,
the taxon_id of the last common ancestor, its scientific/common name and its
lineage.

For example::

    contig_21   172788  uncultured phototrophic eukaryote   cellular organisms,environmental samples

If the *-r* is used, by passing the fasta file containing the nucleotide
sequences the output file is a GFF where for each an annotation for the full
contig length contains the same information of the tab separated file format.

The **lca_line** command accept as input a file where each line consist of a
list of taxon_ids. The separator for the list can be changed and it defaults to
TAB. The last common ancestor for all taxa on a line is searched. The ouput of
this command is the same as the tab separated file of the **lca** command, with
the difference that instead of the first column, which in this command becames
a list of all *taxon_ids* that were used to find the last common ancestor for
that line. The list of *taxon_ids* is separated by semicolon ";".

"""
from __future__ import division
import sys
import argparse
import logging

import mgkit
from . import utils
from mgkit.io import gff, fasta
from mgkit import taxon

LOG = logging.getLogger(__name__)


def set_common_options(parser):
    parser.add_argument(
        '-t',
        '--taxonomy',
        type=str,
        default=None,
        help='Taxonomy file',
        required=True
    )
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


def set_lca_contig_parser(parser):
    parser.add_argument(
        '-b',
        '--bitscore',
        default=0,
        type=float,
        help='Minimum bitscore accepted'
    )
    parser.add_argument(
        '-r',
        '--reference',
        default=None,
        type=argparse.FileType('r'),
        help='Reference file for the GFF, if supplied a GFF file is the output'
    )

    parser.set_defaults(func=lca_contig_command)


def write_lca_gff(file_handle, seq_id, seq, taxon_id, taxon_name, lineage):
    annotation = gff.from_sequence(
        seq_id,
        seq,
        feat_type='LCA',
        taxon_id=taxon_id,
        taxon_name=taxon_name,
        lineage=lineage
    )
    annotation.to_file(file_handle)


def write_lca_tab(file_handle, seq_id, taxon_id, taxon_name, rank, lineage):
    file_handle.write(
        "{}\t{}\t{}\t{}\t{}\n".format(
            seq_id,
            taxon_id,
            taxon_name,
            rank,
            lineage
        )
    )


def get_taxon_info(taxonomy, taxon_id):
    if taxonomy[taxon_id].s_name:
        taxon_name = taxonomy[taxon_id].s_name
    else:
        taxon_name = taxonomy[taxon_id].c_name
    lineage = ','.join(
        tx
        for tx in taxon.get_lineage(taxonomy, taxon_id, names=True)
        if tx
    )
    return taxon_name, lineage


def lca_contig_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    if options.reference is not None:
        seqs = dict(fasta.load_fasta(options.reference))
    else:
        seqs = None

    # groups the annotations by sequence, in case they're not sorted
    annotations = gff.group_annotations(
        (
            annotation
            for annotation in gff.parse_gff(options.input_file)
            # only use annotations whose bitscore pass the filter
            # and have a taxon_id
            if (annotation.bitscore >= options.bitscore) and
               (annotation.taxon_id is not None)
        ),
        lambda annotation: annotation.seq_id
    )
    for seq_id, seq_ann in annotations.iteritems():
        try:
            taxon_id = taxon.last_common_ancestor_multiple(
                taxonomy,
                (annotation.taxon_id for annotation in seq_ann)
            )
        except taxon.NoLcaFound as error:
            LOG.error("No LCA found for %s (%s)", seq_id, error)
            continue

        taxon_name, lineage = get_taxon_info(taxonomy, taxon_id)
        if seqs is not None:
            write_lca_gff(
                options.output_file,
                seq_id,
                seqs[seq_id],
                taxon_id,
                taxon_name,
                lineage
            )
        else:
            write_lca_tab(
                options.output_file,
                seq_id,
                taxon_id,
                taxon_name,
                taxonomy[taxon_id].rank,
                lineage
            )


def set_lca_line_parser(parser):
    parser.add_argument(
        '-s',
        '--separator',
        default='\t',
        type=str,
        help='separator for taxon_ids (defaults to TAB)'
    )
    parser.set_defaults(func=lca_line_command)


def lca_line_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    for line in options.input_file:
        taxon_ids = [
            int(taxon_id)
            for taxon_id in line.strip().split(options.separator)
        ]
        try:
            taxon_id = taxon.last_common_ancestor_multiple(
                taxonomy,
                taxon_ids
            )
        except taxon.NoLcaFound as error:
            LOG.error("No LCA found for %s (%s)", taxon_ids, error)
            continue
        taxon_name, lineage = get_taxon_info(taxonomy, taxon_id)
        write_lca_tab(
            options.output_file,
            ';'.join(str(x) for x in taxon_ids),
            taxon_id,
            taxon_name,
            taxonomy[taxon_id].rank,
            lineage
        )


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Taxonomy Utilities',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_lca_contig = subparsers.add_parser(
        'lca',
        help='Finds the last common ancestor for each sequence in a GFF file'
    )

    set_lca_contig_parser(parser_lca_contig)
    set_common_options(parser_lca_contig)
    utils.add_basic_options(parser_lca_contig)

    parser_lca_line = subparsers.add_parser(
        'lca_line',
        help='Finds the last common ancestor for all IDs in a text file line'
    )

    set_lca_line_parser(parser_lca_line)
    set_common_options(parser_lca_line)
    utils.add_basic_options(parser_lca_line)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
