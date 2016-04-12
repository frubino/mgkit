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

.. note::

    Both also accept the *-n* option, to report the config/line and the
    taxon_ids that had no common ancestors. These are treated as errors and do
    not appear in the output file. The no

Filter by Taxon
***************

The **filter** command of this script allows to filter a GFF file using the
*taxon_id* attribute to include only some annotations, or exclude some. The
filter is based on the `mgkit.taxon.is_ancestor` function, and the
`mgkit.filter.taxon.filter_taxon_by_id_list`. It allows to pass a list of
taxon_id (or taxon_names) to the script. The *include* filter will only output
annotations that have one of the passed taxa as ancestors, while the *exclude*
filter will remove those annotations, that have the passed taxa as ancestors,
from the output.

A list of comma separated taxon_ids can be supplied, as for the names. If any
of the the supplied names have multiple taxon_id (e.g. Actinobacteria) the
script exits and in the log can be found the list of duplicates. For cases like
this, it's preferred for the user to supply a taxon_id, as they can be searched
in NCBI taxonomy (also Uniprot).

.. warning::

    Annotations with no taxon_id are not included in the output of both filters

"""
from __future__ import division
import sys
import argparse
import logging
import functools

import mgkit
from . import utils
from mgkit.io import gff, fasta
from mgkit import taxon
from mgkit.simple_cache import memoize
from mgkit.filter.taxon import filter_taxon_by_id_list

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


def lca_options(parser):
    parser.add_argument(
        '-n',
        '--no-lca',
        type=argparse.FileType('w'),
        default=None,
        help='File to which write records with no LCA',
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
    parser.add_argument(
        '-s',
        '--sorted',
        action='store_true',
        help='''If the GFF file is sorted (all of a sequence annotations are
                contiguos) can use less memory, `sort -s -k 1,1` can be
                used''',
        default=False
    )

    lca_options(parser)

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
        for tx in taxon.get_lineage(taxonomy, taxon_id, names=True, only_ranked=True)
        if tx
    )
    return taxon_name, lineage


def write_no_lca(file_handle, seq_id, taxon_ids):
    file_handle.write(
        "{}\t{}\n".format(
            seq_id,
            ','.join(str(taxon_id) for taxon_id in set(taxon_ids))
        )
    )


def lca_contig_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    if options.no_lca is not None:
        LOG.info(
            "Writing contigs/taxon_ids of contigs with no LCA to (%s)",
            getattr(options.no_lca, 'name', repr(options.no_lca))
        )

    taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    if options.reference is not None:
        seqs = dict(fasta.load_fasta(options.reference))
    else:
        seqs = None

    # basic filter for the presence of a taxon_id and bitscore
    annotations_iter = (
        annotation
        for annotation in gff.parse_gff(options.input_file)
        # only use annotations whose bitscore pass the filter
        # and have a taxon_id
        if (annotation.bitscore >= options.bitscore) and
           (annotation.taxon_id is not None)
    )

    if options.sorted:
        LOG.info("Input GFF is assumed sorted")
        annotations = gff.group_annotations_sorted(
            annotations_iter,
            lambda annotation: annotation.seq_id
        )
    else:
        # groups the annotations by sequence, in case they're not sorted
        annotations = gff.group_annotations(
            annotations_iter,
            lambda annotation: annotation.seq_id
        ).itervalues()

    for seq_ann in annotations:
        seq_id = seq_ann[0].seq_id
        try:
            taxon_id = taxon.last_common_ancestor_multiple(
                taxonomy,
                (annotation.taxon_id for annotation in seq_ann)
            )
        except taxon.NoLcaFound as error:
            LOG.error("No LCA found for %s (%s)", seq_id, error)
            if options.no_lca is not None:
                write_no_lca(
                    options.no_lca,
                    seq_id,
                    (annotation.taxon_id for annotation in seq_ann)
                )
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
    lca_options(parser)
    parser.set_defaults(func=lca_line_command)


def lca_line_command(options):
    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    if options.no_lca is not None:
        LOG.info(
            "Writing taxon_ids of lines with no LCA to (%s)",
            getattr(options.no_lca, 'name', repr(options.no_lca))
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
            if options.no_lca is not None:
                write_no_lca(
                    options.no_lca,
                    '',
                    taxon_ids
                )
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


def validate_taxon_ids(taxon_ids, taxonomy):
    taxon_ids = set(
        int(taxon_id)
        for taxon_id in taxon_ids.split(',')
        if taxon_id
    )
    broken = False
    for taxon_id in taxon_ids:
        if taxon_id not in taxonomy:
            LOG.critical(
                "taxon_id %d was not found in the taxonomy",
                taxon_id
            )
            broken = True
    if broken:
        utils.exit_script(
            "Some of taxon_ids were not found", 1
        )
    return taxon_ids


def validate_taxon_names(taxon_names, taxonomy):
    taxon_names = tuple(
        taxon_name
        for taxon_name in taxon_names.split(',')
        if taxon_name
    )

    broken = False
    taxon_ids = set()
    for taxon_name in taxon_names:
        try:
            taxon_id = taxonomy.find_by_name(taxon_name)
        except KeyError:
            LOG.critical("'%s' was not found in the taxonomy", taxon_name)
            broken = True
            continue
        if len(taxon_id) > 1:
            LOG.critical(
                "'%s' was found multiple times in the taxonomy: (%s)",
                taxon_name,
                ', '.join(str(tid) for tid in taxon_id)
            )
            broken = True
            continue
        # takes the first (and only) taxon_id from the list returned
        taxon_ids.add(taxon_id[0])

    if broken:
        utils.exit_script(
            "Some of taxon_names were not found", 2
        )
    return taxon_ids


def set_filter_taxa_parser(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-i',
        '--include-taxon-id',
        action='store',
        type=str,
        help='Include only taxon_ids (comma separated)'
    )
    group.add_argument(
        '-in',
        '--include-taxon-name',
        action='store',
        type=str,
        help='Include only taxon_names (comma separated)'
    )
    group.add_argument(
        '-e',
        '--exclude-taxon-id',
        action='store',
        type=str,
        help='Exclude taxon_ids (comma separated)'
    )
    group.add_argument(
        '-en',
        '--exclude-taxon-name',
        action='store',
        type=str,
        help='Exclude taxon_names (comma separated)'
    )
    parser.set_defaults(func=filter_taxa_command)


def filter_taxa_command(options):
    taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    exclude = False

    if options.exclude_taxon_name is not None:
        exclude = True
        taxon_ids = validate_taxon_names(
            options.exclude_taxon_name,
            taxonomy
        )
    elif options.exclude_taxon_id is not None:
        exclude = True
        taxon_ids = validate_taxon_ids(
            options.exclude_taxon_id,
            taxonomy
        )
    elif options.include_taxon_name is not None:
        taxon_ids = validate_taxon_names(
            options.include_taxon_name,
            taxonomy
        )
    elif options.include_taxon_id is not None:
        taxon_ids = validate_taxon_ids(
            options.include_taxon_id,
            taxonomy
        )

    if exclude:
        LOG.info("Excluding Taxa")
    else:
        LOG.info("Only include Taxa")

    comp_taxa = functools.partial(
        filter_taxon_by_id_list,
        filter_list=taxon_ids,
        exclude=exclude,
        func=functools.partial(
            taxon.is_ancestor,
            taxonomy
        )
    )
    comp_taxa = memoize(comp_taxa)

    for annotation in gff.parse_gff(options.input_file):
        if annotation.taxon_id is None:
            continue
        if comp_taxa(annotation.taxon_id):
            annotation.to_file(options.output_file)


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

    parser_filter_taxa = subparsers.add_parser(
        'filter',
        help='Filter a GFF file based on taxonomy'
    )

    set_filter_taxa_parser(parser_filter_taxa)
    set_common_options(parser_filter_taxa)
    utils.add_basic_options(parser_filter_taxa)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    mgkit.logger.config_log(options.verbose)
    options.func(options)
