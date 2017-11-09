"""
Add more information to GFF annotations: gene mappings, coverage, taxonomy,
etc..

Uniprot Command
***************

If the *gene_id* of an annotation is a Uniprot ID, the script queries Uniprot
for the requested information. At the moment the information that can be added
is the taxon_id, taxon_name, lineage and mapping to EC, KO, eggNOG IDs.

It's also possible to add mappings to other databases using the *-m* option
with the correct identifier for the mapping, which can be found at `this page
<http://www.uniprot.org/faq/28>`_; for example if it's we want to add the
mappings of uniprot IDs to *BioCyc*, in the *abbreviation* column of the
mappings we find that it's identifier is *REACTOME_ID*, so we pass
*-m REACTOME* to the script (leaving *_ID* out). Mapped IDs are separated by
commas.

The taxonomy IDs are not overwritten if they are found in the annotations, the
*-f* is provided to force the overwriting of those values.

See also :ref:`gff-specs` for more informations about the GFF specifications
used.

.. note::

    As the script needs to query Uniprot a lot, it is recommended to split
    the GFF in several files, so an error in the connection doesn't waste time.

    However, a cache is kept to reduce the number of connections

Taxonomy Command
****************

To refine the taxonomic assignments of predicted genes annotations, the
annotation sequences may be searched against a database like the NCBI *nt*.

This commands takes as input a GFF file, one or more blast output files and a
file with all mappings from GIDs to taxonomy IDs. More information on how to
get the file can be read in the documentation of the function
:func:`mgkit.io.blast.parse_gi_taxa_table`.

The fasta sequences used with BLAST must have as name the uid of the
annotations they refer to, and one way to obtain these sequences is to use the
function :func:`mgkit.io.gff.extract_nuc_seqs` and save them to a fasta file.
Another options is to use the `sequence` command of the `get-gff-info` script
(:ref:`get-gff-info`).

The command accept a minimum bitscore to accept an hit and the taxon ID is
selected by default using top hit method, but LCA can be used, using the *-l*
switch.

Top Hit
+++++++

The best hit is selected from all those found for a sequence which has the
maximum bitscore and identity, with the bitscore having the highest priority.

LCA Taxon
+++++++++

Activated with the *-l* switch, it selects the last common ancestor of all
taxon IDs that are from the cellular organism root in the taxonomy and are
within a 10 bits (by default, can be customised with *-a*) from the hit with
the highest bitscore. If a taxon ID is not found in the taxonomy, it is
excluded. One of the requirements of this option is a file that contains the
full taxonomy from Uniprot/NCBI. The file can be obtained with the following
command::

    $ download_data -x -p -m your@email

The command will output a `taxonomy.pickle` file that can be passed to the `-x`
option :ref:`download-data`.

Coverage Command
****************

Adds coverage information from BAM alignment files to a GFF file, using the
function :func:`mgkit.align.add_coverage_info`, the user needs to supply for
each sample a BAM file, using the `-a` option, whose parameter is in the form
`sample,samplealg.bam`. More samples can be supplied adding more `-a`
arguments.

.. hint::

    As an example, to add coverage for `sample1`, `sample2` the command line
    is::

        add-gff-info coverage -a sample1,sample1.bam -a sample2,sample2.bam \\
        inputgff outputgff

A total coverage for the annotation is also calculated and stored in the
`cov` attribute, while each sample coverage is stored into `sample_cov` as per
:ref:`gff-specs`.

Adding Coverage from samtools depth
***********************************

The *cov_samtools* allows the use of the output of *samtools* **depth**
command. The *-aa* options must be used to pass information about all base
pairs and sequences coverage in the BAM/SAM file. The command accept only one
sample and the relative file/stream from *samtools*, meaning that multiple
samples coverage information must be added one at a time. One solution is to
pipe multiple commands to obtain result wanted. For example::

    $ add-gff-info cov_samtools -s SAMPLE1 -d sample1-coverage input.gff | \
        add-gff-info cov_samtools -s SAMPLE2 -d sample2-coverage - output.gff

This command will add the coverage information for SAMPLE1 and SAMPLE2 from
the respective files.

Uniprot Offline Mappings
************************

Similar to the *uniprot* command, it uses the `idmapping <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_
file provided by Uniprot, which speeds up the process of adding mappings and
taxonomy IDs from Uniprot gene IDs. It's not possible tough to add *EC*
mappings with this command, as those are not included in the file.

Kegg Information
****************

The *kegg* command allows to add information to each annotation. Right now the
information that can be added is restricted to the pathway(s) (reference KO) a
KO is part of and both the KO and pathway(s) descriptions. This information is
stored in keys starting with **ko_**.

Expected Aminoacidic Changes
****************************

Some scripts, like :ref:`snp-parser`, require information about the expected
number of synonymous and non-synonymous changes of an annotation. This can be
done using :meth:`mgkit.io.gff.Annotation.add_exp_syn_count` by the user of the
command `exp_syn` of this script. The attributes added to each annotation are
explained in the :ref:`gff-specs`

Adding Information from eggNOG
******************************

The *eggnog* command allows to add information from the *annotations* file
available for profiles in eggNOG.

Adding Count Data
*****************

Count data on a per-sample basis can be added with the *counts* command. The
accepted inputs are from HTSeq-count and featureCounts. The ouput produced by
featureCounts, is the one from using its **-f** option must be used.

This script accept by default a tab separated file, with a uid in the first
column and the other columns are the counts for each sample, in the same order
as they are passed to the **-s** option. To use the featureCounts file format,
this script **-e** option must be used.

The sample names must be provided in the same order as the columns in the input
files. If the counts are FPKMS the *-f* option can be used.

Adding Taxonomy from a Table
****************************

There are cases where it may needed or preferred to add the taxonomy from a
*gene_id* already provided in the GFF file. For such cases the *addtaxa*
command can be used. It works in a similar way to the *taxonomy* command, only
it expect three different type of inputs:

    * *GI-Taxa* table from NCBI (e.g. gi_taxid_nucl.dmp, )
    * tab separated table
    * dictionary
    * HDF5

The first two are tab separated files, where on each line, the first column is
the *gene_id* that is found in the first column, while the second if the
*taxon_id*.

The third option is a serialised Python *dict*/hash table, whose keys are the
*gene_id* and the value is that gene corresponding *taxon_id*. The serialised
formats accepted are msgpack, json and pickle. The *msgpack* module must be
importable. The option to use json and msgpack allow to integrate this script
with other languages without resorting to a text file.

The last option is a HDF5 created using the *to_hdf* command in
:ref:`taxon-utils`. This requires `pandas` installed and `pytables` and it
provides faster lookup of IDs in the table.

While the default is to look for the *gene_id* attribute in the GFF annotation,
another attribute can be specified, using the **-gene-attr** option.

.. note::

    the dictionary content is loaded after the table files and its keys and
    corresponding values takes precedence over the text files.

.. warning::

    from September 2016 NCBI will retire the GI. In that case the same
    kind of table can be built from the *nucl_gb.accession2taxid.gz* file
    The format is different, but some information can be found in
    :func:`mgkit.io.blast.parse_accession_taxa_table`


Adding information from Pfam
****************************

Adds the Pfam description for the annotation, by downloading the list from
Pfam.

The options allow to specify in which attribute the ID/ACCESSION is stored
(defaults to *gene_id*) and which one between ID/ACCESSION is the value of that
attribute (defaults to *ID*). if no description is found for the family, a
warning message is logged.

Changes
*******

.. versionchanged:: 0.3.3
    changed how *addtaxa* *-a* works, to allow the use of *seq_id* as key to
    add the taxon_id

.. versionchanged:: 0.3.0
    added *cov_samtools* command, *--split* option to *exp_syn*, *-c* option to
    *addtaxa*

.. versionchanged:: 0.2.6
    added *skip-no-taxon* option to *addtaxa*

.. versionchanged:: 0.2.5
    if a dictionary is supplied to *addtaxa*, the GFF is not preloaded

.. versionchanged:: 0.2.3
    added *pfam* command, renamed *gitaxa* to *addtaxa* and made it general

.. versionchanged:: 0.2.2
    added *eggnog*, *gitaxa* and *counts* command

.. versionchanged:: 0.2.1

* added *-d* to *uniprot* command
* added cache to *uniprot* command
* added *kegg* command (cached)

.. versionchanged:: 0.1.16
    added *exp_syn* command

.. versionchanged:: 0.1.15
    *taxonomy* command *-b* option changed

.. versionchanged:: 0.1.13

* added *--force-taxon-id* option to the *uniprot* command
* added *coverage* command
* added *taxonomy* command
* added *unipfile* command

.. versionadded:: 0.1.12
"""
from __future__ import division
import sys
import argparse
import logging
import itertools
import functools
import pysam
import json
import pickle
import progressbar
import mgkit
from . import utils
from .. import align
from .. import logger
from .. import taxon
from .. import kegg
from .. import DependencyError
from ..io import gff, blast, fasta, compressed_handle, open_file
from ..io import uniprot as uniprot_io
from ..net import uniprot as uniprot_net
from ..net import pfam
from ..utils.dictionary import cache_dict_file, HDFDict


LOG = logging.getLogger(__name__)


def set_common_options(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input GFF file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Output GFF file, defaults to stdout'
    )


def set_kegg_parser(parser):
    parser.add_argument(
        '-c',
        '--email',
        action='store',
        type=str,
        help='Contact email',
        default=None
    )
    group = parser.add_argument_group('Requires Internet connection')
    group.add_argument(
        '-d',
        '--description',
        action='store_true',
        default=False,
        help='Add Kegg description'
    )
    group.add_argument(
        '-p',
        '--pathways',
        action='store_true',
        default=False,
        help='Add pathways ID involved'
    )
    group.add_argument(
        '-m',
        '--kegg-id',
        action='store',
        default='gene_id',
        help="""
        In which attribute the Kegg ID is stored (defaults to *gene_id*)"""
    )
    parser.set_defaults(func=kegg_command)


def kegg_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    kegg_client = kegg.KeggClientRest()

    LOG.info('Retrieving KO names')
    ko_names = kegg_client.get_names('ko')

    if options.pathways:
        LOG.info('Retrieving Pathways names')
        # Changes the names of the keys to *ko* instead of *map*
        path_names = dict(
            (path_id.replace('map', 'ko'), name)
            for path_id, name in kegg_client.get_names('path').iteritems()
        )

    ko_cache = {}

    for annotation in gff.parse_gff(options.input_file, gff_type=gff.from_gff):
        try:
            ko_id = annotation.attr[options.kegg_id]
        except KeyError:
            continue
        # if more than one KO is defined
        if ',' in ko_id:
            LOG.warning(
                'More than one KO assigned, skipping annotation: %s',
                annotation.uid
            )
            annotation.to_file(options.output_file)
            continue
        try:
            ko_info = ko_cache[ko_id]
        except KeyError:
            ko_info = {}
            if options.pathways:
                ko_info['pathway'] = kegg_client.link_ids('path', ko_id)
                # If left empty, no pathway will be saved
                if ko_info['pathway']:
                    ko_info['pathway'] = list(
                        path_id
                        for path_id in ko_info['pathway'][ko_id]
                        if path_id.startswith('ko')
                    )
            ko_cache[ko_id] = ko_info
        if options.description:
            annotation.attr['ko_description'] = ko_names.get(ko_id, ko_id)
        if options.pathways and ko_info['pathway']:
            annotation.attr['ko_pathway'] = ','.join(ko_info['pathway'])
            annotation.attr['ko_pathway_names'] = ','.join(
                path_names[path_id]
                for path_id in ko_info['pathway']
            )

        annotation.to_file(options.output_file)


def set_uniprot_parser(parser):
    parser.add_argument(
        '-c',
        '--email',
        action='store',
        type=str,
        help='Contact email',
        default=None
    )
    parser.add_argument(
        '--buffer',
        action='store',
        type=int,
        help='Number of annotations to keep in memory',
        default=50
    )
    group = parser.add_argument_group('Requires Internet connection')
    group.add_argument(
        '-f',
        '--force-taxon-id',
        action='store_true',
        default=False,
        help='Overwrite taxon_id if already present'
    )
    group.add_argument(
        '-t',
        '--taxon-id',
        action='store_true',
        default=False,
        help="""
             Add taxonomic ids to annotations, if taxon_id is found, it won't
             be Overwritten.
             """
    )
    group.add_argument(
        '-l',
        '--lineage',
        action='store_true',
        default=False,
        help='Add taxonomic lineage to annotations'
    )
    group.add_argument(
        '-e',
        '--eggnog',
        action='store_true',
        default=False,
        help='Add eggNOG mappings to annotations'
    )
    group.add_argument(
        '-ec',
        action='store_true',
        default=False,
        help='Add EC mappings to annotations'
    )
    group.add_argument(
        '-ko',
        action='store_true',
        default=False,
        help='Add KO mappings to annotations'
    )
    group.add_argument(
        '-d',
        '--protein-names',
        action='store_true',
        default=False,
        help='Add Uniprot description'
    )
    group.add_argument(
        '-m',
        '--mapping',
        action='append',
        type=str,
        help='Add any DB mappings to annotations'
    )

    parser.set_defaults(func=uniprot_command)


def add_uniprot_info(annotations, options, info_cache):
    columns = []

    if options.taxon_id:
        columns.append('organism')
        columns.append('organism-id')
    if options.eggnog:
        columns.append('database(EGGNOG)')
    if options.ko:
        columns.append('database(KO)')
    if options.ec:
        columns.append('ec')
    if options.lineage:
        columns.append('lineage(ALL)')
    if options.mapping is not None:
        for db in options.mapping:
            columns.append('database({0})'.format(db))
    if options.protein_names:
        columns.append('protein names')

    if not columns:
        return

    LOG.info("Retrieving gene information from Uniprot")

    gene_ids = set(
        x.gene_id
        for x in annotations
        if x.gene_id not in info_cache
    )

    # It may be that all gene_ids are already cached. In which case there's no
    # point querying Uniprot (it will probably give an error)
    if gene_ids:
        data = uniprot_net.get_gene_info(
            gene_ids,
            columns=columns,
            contact=options.email
        )

        info_cache.update(data)

    for annotation in annotations:
        try:
            gene_info = info_cache[annotation.gene_id]
        except KeyError:
            # no data was found
            continue

        for column, values in gene_info.iteritems():
            # nothing found
            if not values:
                continue

            if column == 'organism-id':
                if (annotation.taxon_id and options.force_taxon_id) or \
                   (annotation.taxon_id is None):
                    annotation.attr['taxon_id'] = int(values)
                    annotation.attr['taxon_db'] = 'UNIPROT'
                    # test with a try/expect maybe
                    if 'organism' in columns:
                        annotation.attr['taxon_name'] = gene_info['organism']
            elif column.startswith('lineage'):
                    annotation.attr['lineage'] = gene_info['lineage(ALL)']
            elif column.startswith('database'):
                annotation.attr[
                    'map_{0}'.format(
                        column[:-1].split('(')[1]
                    )
                ] = ','.join(values)
            elif column == 'ec':
                if isinstance(values, list):
                    annotation.attr['EC'] = ','.join(
                        x.split('-')[0]
                        for x in values
                    )
                else:
                    annotation.attr['EC'] = values
            elif column == 'protein names':
                annotation.attr['uniprot_description'] = values


def uniprot_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    if options.buffer < 1:
        options.buffer = 1

    ann_buffer = []

    info_cache = {}

    for annotation in gff.parse_gff(options.input_file, gff_type=gff.from_gff):

        ann_buffer.append(annotation)

        if len(ann_buffer) == options.buffer:

            add_uniprot_info(ann_buffer, options, info_cache)

            for annotation in ann_buffer:
                annotation.to_file(options.output_file)

            ann_buffer = []
    else:
        add_uniprot_info(ann_buffer, options, info_cache)

        for annotation in ann_buffer:
            annotation.to_file(options.output_file)


def set_blast_taxonomy_parser(parser):
    parser.add_argument(
        '-t',
        '--gi-taxa-table',
        action='store',
        default=None,
        required=True,
        help="GIDs taxonomy table (e.g. gi_taxid_nucl.dmp.gz)"
    )
    parser.add_argument(
        '-b',
        '--blast-output',
        action='append',
        default=None,
        required=True,
        help="BLAST output file(s)"
    )
    parser.add_argument(
        '-s',
        '--bitscore',
        action='store',
        default=40,
        type=float,
        help="Minimum bitscore allowed"
    )
    parser.add_argument(
        '-d',
        '--taxon-db',
        action='store',
        default='NCBI-NT',
        help="NCBI database used"
    )
    group = parser.add_argument_group('LCA options')
    group.add_argument(
        '-l',
        '--lca',
        action='store_true',
        default=False,
        help="Use last common ancestor to solve ambiguities"
    )
    group.add_argument(
        '-x',
        '--taxonomy',
        action='store',
        type=argparse.FileType('r'),
        default=None,
        help="Taxonomy file"
    )
    group.add_argument(
        '-a',
        '--max-diff',
        action='store',
        default=10,
        type=float,
        help="Bitscore difference from the max hit"
    )
    parser.set_defaults(func=taxonomy_command)


def get_gids(uid_gid_map):

    gids = set()

    for hits in uid_gid_map.itervalues():
        gids.update(x[0] for x in hits)

    return gids


def choose_by_lca(hits, taxonomy, gid_taxon_map, score=10):
    # the minimum score required to be part of the taxon IDs selected is to
    # be at most 10 bits (by default) from the maximum bitscore found in the
    # hits
    min_score = max(hits, key=lambda x: x[-1])[-1] - score

    # filter over only cellular organisms
    cellular_organism = 131567

    taxon_ids = set()

    for gid, identity, bitscore in hits:
        # No taxon_id found for the gid
        taxon_id = gid_taxon_map.get(gid, None)

        # skip already included taxon_ids
        if taxon_id in taxon_ids:
            continue

        # if the requirements are met
        if (taxon_id is not None) and (bitscore >= min_score):
            try:
                # test to make sure it is a cellular organism
                if taxon.is_ancestor(taxonomy, taxon_id, cellular_organism):
                    taxon_ids.add(taxon_id)
                # not part of the loop
                else:
                    LOG.debug("%d not part of cellular organisms", taxon_id)
            # The taxon id is not in the taxonomy used. It's skipped
            except KeyError:
                LOG.warning("%d is not part of the taxonomy", taxon_id)
                continue

    # no taxon_id passes the filters
    if len(taxon_ids) == 0:
        return None

    # log message fired up only if the package is set on DEBUG
    if mgkit.DEBUG and (len(taxon_ids) > 1):
        LOG.debug(
            "Num hits %d: %s",
            len(taxon_ids),
            ','.join(
                taxonomy[taxon_id].s_name
                for taxon_id in taxon_ids
                if taxon_id in taxonomy
            )
        )

    func = functools.partial(
        taxon.last_common_ancestor,
        taxonomy
    )

    return reduce(func, taxon_ids)


def choose_by_score(hits, gid_taxon_map):
    return gid_taxon_map.get(max(hits, key=lambda x: (x[-1], x[-2]))[0], None)


def taxonomy_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    if options.lca:
        LOG.info(
            "Using LCA to resolve multiple hits %.2f bits from the top hit",
            options.max_diff
        )
        if options.taxonomy is None:
            utils.exit_script('A taxonomy file is required', 1)
        options.taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    uid_gid_map = dict(
        itertools.chain(
            *(
                blast.parse_fragment_blast(x, bitscore=options.bitscore)
                for x in options.blast_output
            )
        )
    )

    gids = get_gids(uid_gid_map)

    gid_taxon_map = dict(
        (gid, taxon_id)
        for gid, taxon_id in blast.parse_gi_taxa_table(
            options.gi_taxa_table, gids=gids
        )
    )

    count = 0
    tot_count = 0

    for annotation in gff.parse_gff(options.input_file):
        tot_count += 1

        hits = uid_gid_map.get(annotation.uid, None)

        if hits is not None:
            if options.lca:
                taxon_id = choose_by_lca(
                    hits,
                    options.taxonomy,
                    gid_taxon_map,
                    options.max_diff
                )
                if taxon_id is not None:
                    LOG.debug(
                        "Selected %s (%s) by LCA",
                        options.taxonomy[taxon_id].s_name,
                        options.taxonomy[taxon_id].rank
                    )
            else:
                taxon_id = choose_by_score(
                    hits,
                    gid_taxon_map
                )
        else:
            taxon_id = None

        if taxon_id is not None:
            count += 1
            annotation.taxon_id = taxon_id
            annotation.taxon_db = options.taxon_db

        annotation.to_file(options.output_file)

    LOG.info(
        "Added taxonomy information to %.2f%% annotations (%d/%d)",
        count / tot_count * 100,
        count,
        tot_count
    )


def split_sample_alg(argument):
    "Split sample/alignment option"
    try:
        sample, bam_file_name = argument.split(',', 1)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Can't get get both sample and bam file from '%s'" % argument
        )

    return (sample, bam_file_name)


def set_coverage_parser(parser):
    parser.add_argument(
        '-a',
        '--sample-alignment',
        action='append',
        required=True,
        type=split_sample_alg,
        help='sample name and correspondent alignment file separated by comma'
    )
    parser.set_defaults(func=coverage_command)


def coverage_command(options):
    samples = []
    bam_files = []

    for sample, bam_file_name in options.sample_alignment:
        samples.append(sample)
        bam_files.append(pysam.Samfile(bam_file_name, 'rb'))

    if len(samples) != len(set(samples)):
        LOG.critical("There are duplicate sample names")
        return 1

    annotations = list(gff.parse_gff(options.input_file))

    align.add_coverage_info(annotations, bam_files, samples)

    gff.write_gff(annotations, options.output_file)


def set_exp_syn_parser(parser):
    """
    .. versionadded:: 0.1.16
    """
    parser.add_argument(
        '-r',
        '--reference',
        required=True,
        type=argparse.FileType('r'),
        help='reference sequence in fasta format'
    )
    parser.add_argument(
        '-s',
        '--split',
        action='store_true',
        help='''Split the sequence header of the reference at the first
        space, to emulate BLAST behaviour''',
        default=False
    )
    parser.set_defaults(func=exp_syn_command)


def exp_syn_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    """
    .. versionadded:: 0.1.16
    """

    seqs = dict(
        (
            seq_id.split(' ')[0] if options.split else seq_id,
            seq
        )
        for seq_id, seq in fasta.load_fasta(options.reference)
    )

    for annotation in gff.parse_gff(options.input_file):
        annotation.add_exp_syn_count(seqs[annotation.seq_id])

        annotation.to_file(options.output_file)


def uniprot_offline_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    LOG.info("Mappings selected: %s", ', '.join(options.mapping))

    annotations = []
    gene_ids = set()

    for annotation in gff.parse_gff(options.input_file):
        annotations.append(annotation)
        gene_ids.add(annotation.gene_id)

    iterator = uniprot_io.uniprot_mappings_to_dict(
        options.mapping_file,
        gene_ids=set(gene_ids),
        mappings=set(options.mapping)
    )

    file_mappings = dict(iterator)

    count = 0

    for annotation in annotations:
        try:
            mappings = file_mappings[annotation.gene_id]
            for mapping_id, mapping_values in mappings.iteritems():
                if mapping_id == uniprot_io.MAPPINGS['taxonomy']:
                    if (annotation.taxon_id and options.force_taxon_id) or \
                       (annotation.taxon_id is None):
                        annotation.taxon_id = mapping_values[0]
                        annotation.taxon_db = 'UNIPROT'
                else:
                    annotation.set_mapping(mapping_id, mapping_values)
            count += 1
        except KeyError:
            pass

        annotation.to_file(options.output_file)

    LOG.info(
        "Number of annotation changed: %d/%d (%.2f%%)",
        count,
        len(annotations),
        count / len(annotations) * 100
    )


def set_uniprot_offline_parser(parser):
    parser.add_argument(
        '-i',
        '--mapping-file',
        action='store',
        default='idmapping.dat.gz',
        required=True,
        help="Uniprot mapping file"
    )
    parser.add_argument(
        '-f',
        '--force-taxon-id',
        action='store_true',
        default=False,
        help="Overwrite taxon_id if already present"
    )
    parser.add_argument(
        '-m',
        '--mapping',
        action='append',
        default=None,
        required=True,
        choices=uniprot_io.MAPPINGS.values(),
        help="Mappings to add"
    )
    parser.set_defaults(func=uniprot_offline_command)


def read_lines_from_files(file_handles):
    file_handles = [
        compressed_handle(file_handle)
        if file_handle.name.endswith('gz')
        else file_handle

        for file_handle in file_handles
    ]
    for line in itertools.chain(*(x for x in file_handles)):
        # the __ is for HTSeq-count files, # for featureCounts
        if line.startswith('#') or line.startswith('__'):
            continue
        yield line.strip()


def eggnog_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    base_info = {}
    LOG.info(
        "Reading eggNOG annotations from files: %s",
        ', '.join(x.name for x in options.annotations_file)
    )
    for line in read_lines_from_files(options.annotations_file):
        level, gene_id, description, source = line.split('\t')
        base_info[gene_id] = dict(
            level=level,
            description=description,
            source=source
        )
    for annotation in gff.parse_gff(options.input_file):
        try:
            ann_info = base_info[annotation.gene_id]
            for key, value in ann_info.iteritems():
                annotation.set_attr(
                    'eggnog_{}'.format(key),
                    value
                )
        except KeyError:
            LOG.warning('No annotation for gene_id %s', annotation.gene_id)
        annotation.to_file(options.output_file)


def set_eggnog_parser(parser):
    parser.add_argument(
        '-a',
        '--annotations-file',
        action='append',
        required=True,
        type=argparse.FileType('r'),
        help="Annotations file"
    )
    parser.set_defaults(func=eggnog_command)


def load_counts(count_files, samples, featureCounts):
    LOG.info("Sample: %s", ', '.join(samples))

    counts = {}

    index = 6 if featureCounts else 1

    for line in read_lines_from_files(count_files):
        # featureCounts puts the header on the second line
        # the first is skipped by read_lines_from_files
        # as it starts with a #
        if featureCounts and line.startswith('Geneid'):
            continue
        line = line.split('\t')
        uid = line[0]
        if len(line[index:]) != len(samples):
            utils.exit_script(
                "The number of samples is not the same as the columns in " +
                "the count file",
                2
            )
        counts[uid] = {}
        for sample, count in zip(samples, line[index:]):
            count = float(count)
            # to save some memory
            if count == 0:
                continue
            counts[uid][sample] = count

    LOG.info("Loaded counts for %s annotations", len(counts))

    return counts


def counts_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    counts = load_counts(
        options.count_files,
        options.samples,
        options.featureCounts
    )

    key = "fpkms_{}" if options.fpkms else "counts_{}"

    for annotation in gff.parse_gff(options.input_file):
        try:
            ann_counts = counts[annotation.uid]
        except KeyError:
            LOG.warning("No counts found for annotation %s", annotation.uid)
            annotation.to_file(options.output_file)
            continue

        for sample in options.samples:
            annotation.attr[key.format(sample)] = ann_counts.get(sample, 0)

        annotation.to_file(options.output_file)


def set_counts_parser(parser):
    parser.add_argument(
        '-s',
        '--samples',
        action='store',
        required=True,
        type=lambda x: [y for y in x.split(',') if y],
        help="Comma separated sample names, in the same order as the count " +
             " file"
    )
    parser.add_argument(
        '-c',
        '--count-files',
        action='append',
        required=True,
        type=argparse.FileType('r'),
        help="Count file(s)"
    )
    parser.add_argument(
        '-f',
        '--fpkms',
        action='store_true',
        default=False,
        help="If the counts are FPKMS"
    )
    parser.add_argument(
        '-e',
        '--featureCounts',
        action='store_true',
        default=False,
        help="If the counts files are from featureCounts (using the -f option)"
    )
    parser.set_defaults(func=counts_command)


def parse_hdf5_arg(argument):
    file_name, table = argument.strip().split(':')
    return (file_name, table)


def addtaxa_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    if (options.gene_taxon_table is not None) and (not options.cache_table):
        annotations = []
        gene_ids = set()

        # the annotations are preloaded if a table is supplied, such as
        # the one from NCBI, so only the gene_ids necessary are taken from that
        # table to save memory

        for annotation in gff.parse_gff(options.input_file):
            annotations.append(annotation)
            gene_ids.add(annotation.attr[options.gene_attr])
        gene_ids = dict(
            blast.parse_accession_taxa_table(
                options.gene_taxon_table,
                acc_ids=gene_ids,
                key=0,
                value=1,
                no_zero=True
            )
        )
    elif (options.gene_taxon_table is not None) and options.cache_table:
        annotations = gff.parse_gff(options.input_file)
        gene_ids = cache_dict_file(
            blast.parse_accession_taxa_table(
                options.gene_taxon_table,
                key=0,
                value=1,
                no_zero=True
            )
        )
    elif options.hdf_table is not None:
        try:
            gene_ids = HDFDict(options.hdf_table[0], options.hdf_table[1])
        except ValueError as e:
            utils.exit_script(
                str(e),
                3
            )
        annotations = gff.parse_gff(options.input_file)
    else:
        # in case a dictionary is supplied, it's expected to fit in memory,
        # meaning that the GFF doesn't have to be preloaded
        annotations = gff.parse_gff(options.input_file)
        gene_ids = {}

    if (options.dictionary is not None) and (options.hdf_table is None):
        dict_file = open_file(options.dictionary, 'r')
        if '.json' in options.dictionary:
            gene_ids.update(json.load(dict_file))
        elif '.msgpack':
            try:
                import msgpack
            except ImportError:
                raise DependencyError('msgpack')
            gene_ids.update(msgpack.load(dict_file))
        else:
            gene_ids.update(pickle.load(dict_file))

        # Ensures all taxon_ids are *int*
        gene_ids = dict(
            (key, int(value))
            for key, value in gene_ids.iteritems()
        )

    if options.taxonomy is not None:
        taxonomy = taxon.UniprotTaxonomy(options.taxonomy)

    if options.hdf_table is not None:
        bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        annotations = bar(annotations)

    for annotation in annotations:
        gene_id = annotation.get_attr(options.gene_attr, str)
        try:
            taxon_id = gene_ids[gene_id]
            annotation.taxon_id = taxon_id
            annotation.set_attr('taxon_db', options.taxon_db)
        except KeyError:
            LOG.error(
                "No Taxon ID for GENE - %s",
                gene_id
            )
            taxon_id = None
            if options.skip_no_taxon:
                continue

        if (options.taxonomy is not None) and (taxon_id is not None):
            try:
                annotation.attr['taxon_name'] = taxonomy[taxon_id].s_name
                annotation.attr['lineage'] = ','.join(
                    taxon_name
                    for taxon_name in taxon.get_lineage(taxonomy, taxon_id, names=True, only_ranked=True)
                    if taxon_name
                )
            except KeyError:
                LOG.warning("Taxon ID %d not found in the Taxonomy", taxon_id)

        annotation.to_file(options.output_file)


def set_addtaxa_parser(parser):
    parser.add_argument(
        '-t',
        '--gene-taxon-table',
        action='store',
        default=None,
        help="""GIDs taxonomy table (e.g. gi_taxid_nucl.dmp.gz) or a similar
                file where GENE/TAXON are tab separated and one per line"""
    )
    parser.add_argument(
        '-f',
        '--hdf-table',
        action='store',
        default=None,
        type=parse_hdf5_arg,
        help="""
        HDF5 file and table name to use for taxon_id lookups. The format to
        pass is the file name, colon and the table file.hf5:taxa-table. The
        index in the table is the accession_id, while the column `taxon_id`
        stores the taxon_id as int
        """
    )
    parser.add_argument(
        '-a',
        '--gene-attr',
        action='store',
        default='gene_id',
        help="""
        In which attribute the GENEID in the table is stored (defaults to
        *gene_id*)"""
    )
    parser.add_argument(
        '-x',
        '--taxonomy',
        action='store',
        default=None,
        help="""
        Taxonomy file - If given, both *taxon_name* and *lineage* attributes
        will be set
        """
    )
    parser.add_argument(
        '-d',
        '--dictionary',
        action='store',
        default=None,
        help="""
        A serialised dictionary, where the key is the GENEID and the value is
        TAXONID. It can be in json or msgpack format (can be a compressed file)
        *Note*: the dictionary values takes precedence over the table files
        """
    )
    parser.add_argument(
        '-e',
        '--skip-no-taxon',
        action='store_true',
        default=False,
        help="""If used, annotations with no taxon_id won't be outputted"""
    )
    parser.add_argument(
        '-db',
        '--taxon-db',
        action='store',
        default='NONE',
        help="DB used to add the taxonomic information"
    )
    parser.add_argument(
        '-c',
        '--cache-table',
        action='store_true',
        default=False,
        help="""If used, annotations are not preloaded, but the taxa table is
        cached, instead.
        """
    )

    parser.set_defaults(func=addtaxa_command)


def pfam_command(options):

    LOG.info(
        'Writing to file (%s)',
        getattr(options.output_file, 'name', repr(options.output_file))
    )

    LOG.info("Downloading Pfam Information")

    pfam_families = pfam.get_pfam_families(
        'acc' if options.use_accession else 'id'
    )

    for annotation in gff.parse_gff(options.input_file):
        pfam_id = annotation.attr[options.id_attr]
        try:
            annotation.attr['pfam_description'] = pfam_families[pfam_id][1]
        except KeyError:
            LOG.warning("No description found for family %s", pfam_id)
        annotation.to_file(options.output_file)


def set_pfam_parser(parser):
    parser.add_argument(
        '-i',
        '--id-attr',
        action='store',
        default='gene_id',
        help="""
        In which attribute the Pfam ID/ACCESSION is stored (defaults to
        *gene_id*)"""
    )
    parser.add_argument(
        '-a',
        '--use-accession',
        action='store_true',
        default=False,
        help="""If used, the attribute value is the Pfam ACCESSION
        (e.g. PF06894), not ID (e.g. Phage_TAC_2)"""
    )
    parser.set_defaults(func=pfam_command)


def set_samtools_depth_parser(parser):
    parser.add_argument(
        '-s',
        '--sample',
        action='store',
        type=str,
        help='sample name'
    )
    parser.add_argument(
        '-d',
        '--depth',
        action='store',
        required=True,
        type=argparse.FileType('r'),
        help='`samtools depth -aa` file'
    )
    parser.add_argument(
        '-n',
        '--num-seqs',
        action='store',
        default=10**4,
        type=int,
        help='Number of sequences to update the log'
    )
    parser.set_defaults(func=samtools_depth_command)


def samtools_depth_command(options):
    depth = align.SamtoolsDepth(options.depth, options.num_seqs)
    if options.sample is None:
        cov_attr = 'cov'
    else:
        cov_attr = '{}_cov'.format(options.sample)

    for annotation in gff.parse_gff(options.input_file):
        cov = depth.region_coverage(
            annotation.seq_id,
            annotation.start,
            annotation.end
        )
        annotation.set_attr(cov_attr, cov)
        annotation.to_file(options.output_file)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Adds informations to a GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    parser_u = subparsers.add_parser(
        'uniprot',
        help='Adds information from GFF whose gene_id is from Uniprot'
    )

    set_uniprot_parser(parser_u)
    set_common_options(parser_u)
    utils.add_basic_options(parser_u, manual=__doc__)

    parser_t = subparsers.add_parser(
        'taxonomy',
        help='''Adds taxonomic information from annotation sequences blasted
                against a NCBI db'''
    )

    set_blast_taxonomy_parser(parser_t)
    set_common_options(parser_t)
    utils.add_basic_options(parser_t, manual=__doc__)

    parser_c = subparsers.add_parser(
        'coverage',
        help='Adds coverage information from BAM Alignment files'
    )

    set_coverage_parser(parser_c)
    set_common_options(parser_c)
    utils.add_basic_options(parser_c, manual=__doc__)

    parser_e = subparsers.add_parser(
        'exp_syn',
        help='Adds expected synonymous and non-synonymous changes information'
    )

    set_exp_syn_parser(parser_e)
    set_common_options(parser_e)
    utils.add_basic_options(parser_e, manual=__doc__)

    parser_f = subparsers.add_parser(
        'unipfile',
        help='Adds mappings and taxonomy from Uniprot mapping file'
    )

    set_uniprot_offline_parser(parser_f)
    set_common_options(parser_f)
    utils.add_basic_options(parser_f, manual=__doc__)

    parser_k = subparsers.add_parser(
        'kegg',
        help='Adds information and mapping from Kegg'
    )

    set_kegg_parser(parser_k)
    set_common_options(parser_k)
    utils.add_basic_options(parser_k, manual=__doc__)

    parser_eggnog = subparsers.add_parser(
        'eggnog',
        help='Adds information from eggNOG'
    )

    set_eggnog_parser(parser_eggnog)
    set_common_options(parser_eggnog)
    utils.add_basic_options(parser_eggnog, manual=__doc__)

    parser_counts = subparsers.add_parser(
        'counts',
        help='Adds counts data to the GFF'
    )

    set_counts_parser(parser_counts)
    set_common_options(parser_counts)
    utils.add_basic_options(parser_counts, manual=__doc__)

    parser_addtaxa = subparsers.add_parser(
        'addtaxa',
        help='''Adds taxonomy information from a GI-Taxa, gene_id/taxon_id
                table or a dictionary serialised as a pickle/msgpack/json file
                '''
    )

    set_addtaxa_parser(parser_addtaxa)
    set_common_options(parser_addtaxa)
    utils.add_basic_options(parser_addtaxa, manual=__doc__)

    parser_pfam = subparsers.add_parser(
        'pfam',
        help='Adds information from Pfam'
    )

    set_pfam_parser(parser_pfam)
    set_common_options(parser_pfam)
    utils.add_basic_options(parser_pfam, manual=__doc__)

    parser_samtools_depth = subparsers.add_parser(
        'cov_samtools',
        help='Adds information from samtools_depth'
    )

    set_samtools_depth_parser(parser_samtools_depth)
    set_common_options(parser_samtools_depth)
    utils.add_basic_options(parser_samtools_depth, manual=__doc__)

    # top parser
    utils.add_basic_options(parser, manual=__doc__)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)
