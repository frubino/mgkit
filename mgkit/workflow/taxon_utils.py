"""
The script contains commands used to access functionality related to
taxonomy, without the need to write ad-hoc code for functionality that
can be part of a workflow. One example is access to the the last common
ancestor function contained in the :mod:`mgkit.taxon`.

Last Common Ancestor (lca and lca_line)
***************************************

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
    not appear in the output file.

Krona Output
############

.. versionadded:: 0.3.0

The *lca* command supports the writing of a file compatible with Krona. The
output file can be used with the *ktImportText/ImportText.pl* script included
with `KronaTools <https://github.com/marbl/Krona/wiki>`_. Specifically, the
output from *taxon_utils* will be a file with all the lineages found (tab
separated), that can be used with::

    $ ktImportText -q taxon_utils_ouput

Note the use of *-q* to make the script count the lineages. Sequences with no
LCA found will be marked as *No LCA* in the graph, the *-n* is not required.

.. note::

    Please note that the output won't include any sequence that didn't have a
    hit with the software used. If that's important, the **-kt** option can be
    used to add a number of *Unknown* lines at the end, to read the total
    supplied.

Filter by Taxon
***************

The **filter** command of this script allows to filter a GFF file using the
*taxon_id* attribute to include only some annotations, or exclude some. The
filter is based on the `mgkit.taxon.is_ancestor` function, and the
`mgkit.filter.taxon.filter_taxon_by_id_list`. It can also filter a table (tab
separated values) when the first element is an ID and the second is a taxon_id.
An example of a table of this sort is the output of the `download-ncbi-taxa.sh`
and `download-uniprot-taxa.sh`, where each accession of a database is associated
to a taxon_id.

Multiple taxon_id can be passed, either for inclusion or exclusion. If both
exclusion and inclusion is used, the first check is on the inclusion and then on
the exclusion. In alternative to passing taxon_id, taxon_names can be passed,
with values such as 'cellular organisms' that needs to be quoted. Example::

    $ taxon-utils filter -i 2 -in archaea -en prevotella -t taxonomy.pickle in.gff out.gff

Which will keep only line that are from Bacteria (taxon_id=2) and exclude those
from the genus *Prevotella*. It will be also include Archaea.

Multiple inclusion and exclusion flags can be put::

    $ taxon-utils filter -i 2 -i 2172 -t taxonomy in.gff out.gff

In particular, the inclusion flag is tested first and then the exclusion is
tested. So a line like this one:

.. code-block:: bash

    printf "TEST\\t838\\nTEST\\t1485" | taxon-utils filter -p -t taxonomy.pickle -i 2 -i 1485 -e 838

Will produce **TEST 1485**, because both Prevotella (838) and Clostridium (1485)
are Bacteria (2) OR Prevotella, but Prevotella must be excluded according to
the exclusion option. This line also illustrate that a tab-separated file, where
the second column contains taxon IDs, can be filtered. In particular it can be
applied to files produced by `download-ncbi-taxa.sh` or
`download-uniprot-taxa.sh` (see :ref:`download-taxonomy`).

.. warning::

    Annotations with no taxon_id are not included in the output of both filters

Convert Taxa Tables to HDF5
***************************

This command is used to convert the taxa tables download from Uniprot and NCBI,
using the scripts mentioned in :ref:`download-data`,
`download-uniprot-taxa.sh` and `download-ncbi-taxa` into a HDF5 file that can
be used with the *addtaxa* command in :ref:`add-gff-info`.

The advantage is a faster lookup of the IDs. The other is a smaller memory
footprint when a great number of annotations are kept in memory.

Changes
*******

.. versionchanged:: 0.3.4
    changed interface and behaviour for *filter*, also now can filter tables;
    *lca* has changed the interface and allows the output of a 2 column table

.. versionchanged:: 0.3.1
    added *to_hdf* command

.. versionchanged:: 0.3.1
    added *-j* option to *lca*, which outputs a JSON file with the LCA results

.. versionchanged:: 0.3.0
    added *-k* and *-kt* options for Krona output, lineage now includes the LCA
    also added *-a* option to select between lineages with only ranked taxa.
    Now it defaults to all components.

.. versionchanged:: 0.2.6
    added *feat-type* option to *lca* command, added phylum output to nolca

.. versionadded:: 0.2.5

"""
from __future__ import division
from builtins import range
from future.utils import viewvalues
import logging
import functools
import json
import click
import mgkit
import pandas as pd
from tqdm import tqdm
from . import utils
from mgkit.io import gff, fasta, blast
from mgkit import taxon
from mgkit.simple_cache import memoize
from mgkit.filter.taxon import filter_taxon_by_id_list


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


def write_lca_gff(file_handle, seq_id, seq, taxon_id, taxon_name, lineage,
                  feat_type):
    annotation = gff.from_sequence(
        seq_id,
        seq,
        feat_type=feat_type,
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
        ).encode('ascii')
    )


def write_lca_tab_simple(file_handle, seq_id, taxon_id):
    file_handle.write(
        "{}\t{}\n".format(seq_id, taxon_id).encode('ascii')
    )


def get_taxon_info(taxonomy, taxon_id, only_ranked):
    if taxonomy[taxon_id].s_name:
        taxon_name = taxonomy[taxon_id].s_name
    else:
        taxon_name = taxonomy[taxon_id].c_name
    lineage = ','.join(
        tx
        for tx in taxon.get_lineage(taxonomy, taxon_id, names=True, only_ranked=only_ranked, with_last=True)
        if tx
    )
    return taxon_name, lineage


def write_no_lca(file_handle, seq_id, taxon_ids, extra=None):
    file_handle.write(
        "{}\t{}\t{}\n".format(
            seq_id,
            ','.join(str(taxon_id) for taxon_id in set(taxon_ids)),
            '' if extra is None else ','.join(extra)
        ).encode('ascii')
    )


def write_krona(file_handle, taxonomy, taxon_id, only_ranked):
    if taxon_id is None:
        lineage = ['No LCA']
    else:
        lineage = taxon.get_lineage(
            taxonomy,
            taxon_id,
            names=True,
            only_ranked=only_ranked,
            with_last=True
        )
    file_handle.write(
        '{}\n'.format('\t'.join(lineage)).encode('ascii')
    )


def write_json(lca_dict, seq_id, taxonomy, taxon_id, only_ranked):
    if taxon_id is None:
        lineage = ""
    else:
        lineage = taxon.get_lineage(
            taxonomy,
            taxon_id,
            names=True,
            only_ranked=only_ranked,
            with_last=True
        )

    lca_dict[seq_id] = dict(
        taxon_id=taxon_id,
        lineage=';'.join(lineage),
    )


@main.command('lca', help="""Finds the last common ancestor for each sequence
              in a GFF file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--taxonomy', required=True, help='Taxonomy file')
@click.option('-n', '--no-lca', type=click.File('wb'), default=None,
              help='File to which write records with no LCA')
@click.option('-a', '--only-ranked', is_flag=True, default=False,
              help='''If set, only taxa that have a rank will be used in the lineageself. This is not advised for lineages such as Viruses, where the top levels have no rank''')
@click.option('-b', '--bitscore', default=0, type=click.FLOAT, show_default=True,
              help='Minimum bitscore accepted')
@click.option('-m', '--rename', is_flag=True, default=False,
              help='Emulates BLAST behaviour for headers (keep left of first space)')
@click.option('-s', '--sorted', is_flag=True, default=False,
              help='''If the GFF file is sorted (all of a sequence annotations are contiguos) can use less memory, `sort -s -k 1,1` can be used''')
@click.option('-ft', '--feat-type', default='LCA', show_default=True,
              help='Feature type used if the output is a GFF (default is *LCA*)')
@click.option('-r', '--reference', default=None, type=click.File('rb'),
              help='Reference file for the GFF, if supplied a GFF file is the output')
@click.option('-p', '--simple-table', is_flag=True, default=False,
              help='Uses a 2 column table format (seq_id taxon_id) TAB separated')
@click.option('-kt', '--krona-total', type=click.INT, default=None,
              help='''Total number of raw sequences (used to output correct percentages in Krona''')
@click.option('-f', '--out-format', default='tab', show_default=True,
              type=click.Choice(['krona', 'json', 'tab', 'gff']),
              help='Format of output file')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('gff-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def lca_contig_command(verbose, taxonomy, no_lca, only_ranked, bitscore,
                       rename, sorted, feat_type, reference, simple_table,
                       krona_total, out_format, progress, gff_file,
                       output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    if (out_format == 'gff') and (reference is None):
        utils.exit_script(
            'The output format *gff* requires a FASTA file, passed with -r', 3)

    if no_lca is not None:
        LOG.info(
            "Writing contigs/taxon_ids of contigs with no LCA to (%s)",
            getattr(no_lca, 'name', repr(no_lca))
        )

    taxonomy = taxon.Taxonomy(taxonomy)

    if reference is not None:
        seqs = dict(
            fasta.load_fasta_rename(reference) if rename else fasta.load_fasta(reference)
        )
    else:
        seqs = None

    # basic filter for the presence of a taxon_id and bitscore
    annotations_iter = (
        annotation
        for annotation in gff.parse_gff(gff_file)
        # only use annotations whose bitscore pass the filter
        # and have a taxon_id
        if (
            (annotation.bitscore >= bitscore) or
            (annotation.bitscore is None)
            ) and (annotation.taxon_id is not None) and
            # redundant probably, but used in cases when a taxon_id was deleted
            # from the taxonomy
            (annotation.taxon_id in taxonomy)
    )

    if sorted:
        LOG.info("Input GFF is assumed sorted")
        annotations = gff.group_annotations_sorted(
            annotations_iter,
            lambda annotation: annotation.seq_id
        )
    else:
        # groups the annotations by sequence, in case they're not sorted
        annotations = viewvalues(gff.group_annotations(
            annotations_iter,
            lambda annotation: annotation.seq_id
        ))

    count = 0
    lca_dict = {}
    if progress:
        annotations = tqdm(annotations)
    for seq_ann in annotations:
        count += 1
        seq_id = seq_ann[0].seq_id
        try:
            taxon_id = taxon.last_common_ancestor_multiple(
                taxonomy,
                (annotation.taxon_id for annotation in seq_ann)
            )
        except taxon.NoLcaFound as error:
            LOG.warning("No LCA found for %s (%s)", seq_id, error)
            if no_lca is not None:
                write_no_lca(
                    no_lca,
                    seq_id,
                    (annotation.taxon_id for annotation in seq_ann),
                    extra=set(
                        taxonomy.get_ranked_taxon(
                            annotation.taxon_id,
                            'phylum'
                        ).s_name
                        for annotation in seq_ann
                    )
                )
            if out_format == 'krona':
                write_krona(output_file, taxonomy, None, False)
            continue

        taxon_name, lineage = get_taxon_info(
            taxonomy,
            taxon_id,
            only_ranked
        )
        if out_format == 'gff':
            write_lca_gff(
                output_file,
                seq_id,
                seqs[seq_id],
                taxon_id,
                taxon_name,
                lineage,
                feat_type
            )
        elif out_format == 'krona':
            write_krona(
                output_file,
                taxonomy,
                taxon_id,
                only_ranked
            )
        elif out_format == 'json':
            write_json(
                lca_dict,
                seq_id,
                taxonomy,
                taxon_id,
                only_ranked
            )
        elif out_format == 'tab':
            if simple_table:
                write_lca_tab_simple(output_file, seq_id, taxon_id)
            else:
                write_lca_tab(
                    output_file,
                    seq_id,
                    taxon_id,
                    taxon_name,
                    taxonomy[taxon_id].rank,
                    lineage
                )
    if (out_format == 'krona') and (krona_total is not None):
        for index in range(count, krona_total):
            output_file.write('Unknown\n')

    if out_format == 'json':
        output_file.write(json.dumps(lca_dict, indent=4).encode('ascii'))


@main.command('lca_line', help="""Finds the last common ancestor for all IDs in
              a text file line""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--taxonomy', required=True, help='Taxonomy file')
@click.option('-n', '--no-lca', type=click.File('wb'), default=None,
              help='File to which write records with no LCA')
@click.option('-a', '--only-ranked', is_flag=True, default=False,
              help='''If set, only taxa that have a rank will be used in the lineageself. This is not advised for lineages such as Viruses, where the top levels have no rank''')
@click.option('-s', '--separator', default='\t',
              help='separator for taxon_ids (defaults to TAB)')
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def lca_line_command(verbose, taxonomy, no_lca, only_ranked, separator,
                     input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    if no_lca is not None:
        LOG.info(
            "Writing taxon_ids of lines with no LCA to (%s)",
            getattr(no_lca, 'name', repr(no_lca))
        )

    taxonomy = taxon.Taxonomy(taxonomy)

    for line in input_file:
        line = line.decode('ascii')
        taxon_ids = [
            int(taxon_id)
            for taxon_id in line.strip().split(separator)
        ]
        try:
            taxon_id = taxon.last_common_ancestor_multiple(
                taxonomy,
                taxon_ids
            )
        except taxon.NoLcaFound as error:
            LOG.error("No LCA found for %s (%s)", taxon_ids, error)
            if no_lca is not None:
                write_no_lca(
                    no_lca,
                    '',
                    taxon_ids
                )
            continue
        taxon_name, lineage = get_taxon_info(taxonomy, taxon_id, only_ranked)
        write_lca_tab(
            output_file,
            ';'.join(str(x) for x in taxon_ids),
            taxon_id,
            taxon_name,
            taxonomy[taxon_id].rank,
            lineage
        )


def validate_taxon_ids(taxon_ids, taxonomy):
    taxon_ids = set(taxon_ids)
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


@main.command('filter', help="Filter a GFF file or a table based on taxonomy")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-p', '--table', is_flag=True, default=False)
@click.option('-t', '--taxonomy', required=True, help='Taxonomy file')
@click.option('-i', '--include-taxon-id', multiple=True, default=None,
              type=click.INT, help='Include only taxon_ids')
@click.option('-in', '--include-taxon-name', multiple=True, default=None,
              help='Include only taxon_names')
@click.option('-e', '--exclude-taxon-id', multiple=True, default=None,
              type=click.INT, help='Exclude taxon_ids')
@click.option('-en', '--exclude-taxon-name', multiple=True, default=None,
              help='Exclude taxon_names')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def filter_taxa_command(verbose, table, taxonomy, include_taxon_id,
                        include_taxon_name, exclude_taxon_id,
                        exclude_taxon_name, progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    taxonomy = taxon.Taxonomy(taxonomy)

    exclude_ids = validate_taxon_ids(exclude_taxon_id, taxonomy) | \
        validate_taxon_names(exclude_taxon_name, taxonomy)

    include_ids = validate_taxon_ids(include_taxon_id, taxonomy) | \
        validate_taxon_names(include_taxon_name, taxonomy)

    if exclude_ids:
        LOG.info("Excluding Taxa: %s", exclude_ids)
        exclude_func = functools.partial(
            filter_taxon_by_id_list,
            filter_list=exclude_ids,
            exclude=True,
            func=functools.partial(
                taxon.is_ancestor,
                taxonomy
            )
        )
        exclude_func = memoize(exclude_func)
    else:
        exclude_func = None
    if include_ids:
        LOG.info("Only include Taxa: %s", include_ids)
        include_func = functools.partial(
            filter_taxon_by_id_list,
            filter_list=include_ids,
            exclude=False,
            func=functools.partial(
                taxon.is_ancestor,
                taxonomy
            )
        )
        include_func = memoize(include_func)
    else:
        include_func = None

    if table:
        iterator = blast.parse_accession_taxa_table(input_file, key=0, value=1,
                                                    num_lines=None)
        if progress:
            iterator = tqdm(iterator)
        for acc_id, taxon_id in iterator:
            if include_func is not None:
                if not include_func(taxon_id):
                    continue
            if exclude_func is not None:
                if not exclude_func(taxon_id):
                    continue
            output_file.write("{}\t{}\n".format(acc_id, taxon_id).encode('ascii'))
    else:
        iterator = gff.parse_gff(input_file)
        if progress:
            iterator = tqdm(iterator)
        for annotation in iterator:
            if annotation.taxon_id is None:
                continue
            if include_func is not None:
                if not include_func(annotation.taxon_id):
                    continue
            if exclude_func is not None:
                if not exclude_func(annotation.taxon_id):
                    continue
            annotation.to_file(output_file)


@main.command('to_hdf', help="""Convert a taxa table to HDF5, with the input as
              tabular format, defaults to stdin. Output file, defaults to
              (taxa-table.hf5)""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-n', '--table-name', default='taxa', show_default=True,
              help='Name of the table/storage to use')
@click.option('-w', '--overwrite', default=False, is_flag=True,
              help='Overwrite the file, instead of appending to it')
@click.option('-s', '--index-size', default=12, type=click.INT,
              show_default=True,
              help='Maximum number of characters for the gene_id')
@click.option('-c', '--chunk-size', default=5000000, type=click.INT,
              show_default=True,
              help='Chunk size to use when reading the input file')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input_file', type=click.File('rb'), default='-')
@click.argument('output_file', default='taxa-table.hf5')
def taxa_table_command(verbose, table_name, overwrite, index_size, chunk_size,
                       progress, input_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    hdf = pd.HDFStore(
        output_file,
        mode='w' if overwrite else 'a'
    )

    iterator = pd.read_table(
        input_file, header=None,
        squeeze=False,
        index_col=0,
        engine='c',
        chunksize=chunk_size,
        names=['taxon_id']
    )
    LOG.info('Reading Taxa Table from file (%s)', input_file.name)
    LOG.info(
        'Writing HDF5 file (%s) Table (%s) - Overwrite: %s',
        output_file,
        table_name,
        overwrite
    )
    if progress:
        iterator = tqdm(iterator)
    for chunk in iterator:
        hdf.append(
            table_name,
            chunk,
            min_itemsize=index_size,
            index=False,
            # data_columns=False
        )

    LOG.info('Creating Indices')
    hdf.create_table_index(table_name, optlevel=9, kind='full')
    hdf.close()

    LOG.info("It's reccomended to compress the file with `ptrepack`")
    LOG.info(
        "e.g. ptrepack --propindexes --complevel 9 --complib blosc %s:/ taxa-table-compressed.hf5:/",
        output_file
    )
