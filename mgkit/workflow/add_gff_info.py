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
command. The command work similarly to *coverage*, accepting compressed *depth*
files as well. If only one *depth* file is passed and no sample is passed, the
attribute in the GFF will be *cov*, otherwise the attribute will be
*sample1_cov*, *sample2_cov*, etc.

.. note::

    From version 0.4.2 the behaviour is different in the following ways:

    * the GFF file is loaded in memory to get information that helps reduce
        the memory size
    * the depth file is scanned and each time a sequence in the GFF file is
        found, coverage data is added, this way memory is freed more often
    * the GFF file is written at the very end, so chaining multiple commands
        (one for each sample) doesn't give any advantage now
    * there's no need to use the `-aa` option of `samtools depth`, since some
        assumptions are made when scanning the fils (see next point)
    * when the depth file is scanned completely, the sequences in the GFF file
        with coverage information are set coverage 0.0, instead of raising an
        error
    * average coverage is only calculated if requested, with the `-m` option

The GFF file is first loaded in memory to get information about the maximum
array size for each sequence, then the script proceeds for each depth file, to
scan it until a) all sequences in the GFF are found or b) until the end of the
file. If b) is the case, all sequences with not found in the depth file are set
to a coverage of 0.0. A warning with the number of sequences lacking coverage
is printed, for diagnostic purposes. If the the mean coverage is requested with
the `-m` option, coverage is then calculated

To create samtools *depth* files, you can use this command::

    $ samtools depth [bam_file] > [depth_file]

    where [bam_file] is the alignment to be used and [depth_file] where the
    file will be stored (by redirection). You could create a gzipped file, as
    the file can be huge (and in my experience compressing reduces by ~10x the
    file size)::

    $ samtools depth [bam_file] | gzip > [depth_file.gz]

.. warning::

    if version < 0.4.2:

    The *-aa* options must be used to pass information about all base
    pairs and sequences coverage in the BAM/SAM file.

    $ samtools depth -aa bam_file

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

.. versionchanged:: 0.4.2
    added *-m* option for *cov_samtools* command, to calculate the average
    coverage for an annotation (*cov* attribute). Fixed loading compressed
    files

.. versionchanged:: 0.3.4
    removed the *taxonomy* command, since a similar result can be obtained with
    *taxon-utils lca* and *add-gff-info addtaxa*. Removed *eggnog* command and
    added option to verbose the logging in *cov_samtools* (now is quiet), also
    changed the interface

.. versionchanged:: 0.3.3
    changed how *addtaxa* *-a* works, to allow the use of *seq_id* as key to
    add the taxon_id

.. versionchanged:: 0.3.0
    added *cov_samtools* command, *--split* option to *exp_syn*, *-c* option to
    *addtaxa*. *kegg* now does not skip annotations when the attribute is not
    found.

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
from builtins import zip
from future.utils import viewitems
import logging
import itertools
import pysam
import json
import pickle
from tqdm import tqdm
import click
import msgpack
import mgkit
import mgkit.counts.func
from . import utils
from .. import align
from .. import logger
from .. import taxon
from .. import kegg
from ..io import gff, blast, fasta, open_file
from ..io import uniprot as uniprot_io
from ..net import uniprot as uniprot_net
from ..net import pfam
from ..utils.dictionary import cache_dict_file, HDFDict


LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


@main.command('kegg', help="""Adds information and mapping from Kegg""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-c', '--email', help='Contact email', required=True)
@click.option('-d', '--description', default=False, is_flag=True,
              help='Add Kegg description')
@click.option('-p', '--pathways', is_flag=True, default=False,
              help='Add pathways ID involved')
@click.option('-m', '--kegg-id', default='gene_id',
              help="""In which attribute the Kegg ID is stored (defaults to *gene_id*)""")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def kegg_command(verbose, email, description, pathways, kegg_id, input_file,
                 output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    kegg_client = kegg.KeggClientRest()

    LOG.info('Retrieving KO names')
    ko_names = kegg_client.get_ids_names('ko')

    if pathways:
        LOG.info('Retrieving Pathways names')
        # Changes the names of the keys to *ko* instead of *map*
        path_names = dict(
            (path_id.replace('map', 'ko'), name)
            for path_id, name in viewitems(kegg_client.get_ids_names('path'))
        )

    ko_cache = {}

    for annotation in gff.parse_gff(input_file, gff_type=gff.from_gff):
        try:
            ko_id = annotation.attr[kegg_id]
        except KeyError:
            annotation.to_file(output_file)
            continue
        # if more than one KO is defined
        if ',' in ko_id:
            LOG.warning(
                'More than one KO assigned, skipping annotation: %s',
                annotation.uid
            )
            annotation.to_file(output_file)
            continue
        try:
            ko_info = ko_cache[ko_id]
        except KeyError:
            ko_info = {}
            if pathways:
                ko_info['pathway'] = kegg_client.link_ids('path', ko_id)
                # If left empty, no pathway will be saved
                if ko_info['pathway']:
                    ko_info['pathway'] = list(
                        path_id
                        for path_id in ko_info['pathway'][ko_id]
                        if path_id.startswith('ko')
                    )
            ko_cache[ko_id] = ko_info
        if description:
            annotation.attr['ko_description'] = ko_names.get(ko_id, ko_id)
        if pathways and ko_info['pathway']:
            annotation.attr['ko_pathway'] = ','.join(ko_info['pathway'])
            annotation.attr['ko_pathway_names'] = ','.join(
                path_names[path_id]
                for path_id in ko_info['pathway']
            )

        annotation.to_file(output_file)


def add_uniprot_info(annotations, email, force_taxon_id, taxon_id, lineage,
                     eggnog, enzymes, kegg_orthologs, protein_names, mapping,
                     info_cache):
    columns = []

    if taxon_id:
        columns.append('organism')
        columns.append('organism-id')
    if eggnog:
        columns.append('database(EGGNOG)')
    if kegg_orthologs:
        columns.append('database(KO)')
    if enzymes:
        columns.append('ec')
    if lineage:
        columns.append('lineage(ALL)')
    if mapping is not None:
        for db in mapping:
            columns.append('database({0})'.format(db))
    if protein_names:
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
            contact=email
        )

        info_cache.update(data)

    for annotation in annotations:
        try:
            gene_info = info_cache[annotation.gene_id]
        except KeyError:
            # no data was found
            continue

        for column, values in viewitems(gene_info):
            # nothing found
            if not values:
                continue

            if column == 'organism-id':
                if (annotation.taxon_id and force_taxon_id) or \
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


@main.command('uniprot', help="""Adds information from GFF whose gene_id is
              from Uniprot""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-c', '--email', help='Contact email', required=True)
@click.option('--buffer', type=click.INT, default=50, show_default=True,
              help='Number of annotations to keep in memory')
@click.option('-f', '--force-taxon-id', default=False, is_flag=True,
              help='Overwrite taxon_id if already present')
@click.option('-t', '--taxon-id', default=False, is_flag=True,
              help="""Add taxonomic ids to annotations, if taxon_id is found, it won't be Overwritten.""")
@click.option('-l', '--lineage', default=False, is_flag=True,
              help='Add taxonomic lineage to annotations')
@click.option('-e', '--eggnog', default=False, is_flag=True,
              help='Add eggNOG mappings to annotations')
@click.option('-ec', '--enzymes', default=False, is_flag=True,
              help='Add EC mappings to annotations')
@click.option('-ko', '--kegg_orthologs', default=False, is_flag=True,
              help='Add KO mappings to annotations')
@click.option('-d', '--protein-names', default=False, is_flag=True,
              help='Add Uniprot description')
@click.option('-m', '--mapping', multiple=True, default=None,
              help='Add any DB mappings to annotations')
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def uniprot_command(verbose, email, buffer, force_taxon_id, taxon_id, lineage,
                    eggnog, enzymes, kegg_orthologs, protein_names, mapping,
                    input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    if buffer < 1:
        buffer = 1

    ann_buffer = []

    info_cache = {}

    for annotation in gff.parse_gff(input_file, gff_type=gff.from_gff):

        ann_buffer.append(annotation)

        if len(ann_buffer) == buffer:

            add_uniprot_info(ann_buffer, email, force_taxon_id, taxon_id,
                             lineage, eggnog, enzymes, kegg_orthologs,
                             protein_names, mapping, info_cache)

            for annotation in ann_buffer:
                annotation.to_file(output_file)

            ann_buffer = []
    else:
        add_uniprot_info(ann_buffer, email, force_taxon_id, taxon_id, lineage,
                         eggnog, enzymes, kegg_orthologs, protein_names,
                         mapping, info_cache)

        for annotation in ann_buffer:
            annotation.to_file(output_file)


def split_sample_alg(ctx, param, values):
    "Split sample/alignment option"

    new_values = []
    for value in values:
        try:
            sample, bam_file_name = value.split(',', 1)
        except ValueError:
            raise click.BadParameter(
                "Can't get get both sample and bam file from '%s'" % value
            )
        new_values.append((sample, bam_file_name))

    return new_values


@main.command('coverage', help="""Adds coverage information from BAM Alignment
              files""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-a', '--sample-alignment', multiple=True, required=True,
              callback=split_sample_alg,
              help='sample name and correspondent alignment file separated by comma')
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def coverage_command(verbose, sample_alignment, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    samples = []
    bam_files = []

    for sample, bam_file_name in sample_alignment:
        samples.append(sample)
        bam_files.append(pysam.Samfile(bam_file_name, 'rb'))

    if len(samples) != len(set(samples)):
        utils.exit_script("There are duplicate sample names", 1)

    annotations = list(gff.parse_gff(input_file))

    align.add_coverage_info(annotations, bam_files, samples)

    gff.write_gff(annotations, output_file)


@main.command('exp_syn', help="""Adds expected synonymous and non-synonymous
              changes information""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-r', '--reference', required=True, type=click.File('rb'),
              help='reference sequence in fasta format')
@click.option('-s', '--split', default=False, is_flag=True,
              help='''Split the sequence header of the reference at the first space, to emulate BLAST behaviour''')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def exp_syn_command(verbose, reference, split, progress, input_file,
                    output_file):
    """
    .. versionadded:: 0.1.16
    """
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    seqs = dict(
        (
            seq_id.split(' ')[0] if split else seq_id,
            seq
        )
        for seq_id, seq in fasta.load_fasta(reference)
    )

    iterator = gff.parse_gff(input_file)

    if progress:
        iterator = tqdm(iterator)

    for annotation in iterator:
        annotation.add_exp_syn_count(seqs[annotation.seq_id])

        annotation.to_file(output_file)


@main.command('unipfile', help="""Adds expected synonymous and non-synonymous
              changes information""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-i', '--mapping-file', default='idmapping.dat.gz', required=True,
              help="Uniprot mapping file")
@click.option('-f', '--force-taxon-id', default=False, is_flag=True,
              help="Overwrite taxon_id if already present")
@click.option('-m', '--mapping', multiple=True, default=None, required=True,
              type=click.Choice(list(uniprot_io.MAPPINGS.values())),
              help="Mappings to add")
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def uniprot_offline_command(verbose, mapping_file, force_taxon_id, mapping,
                            progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    LOG.info("Mappings selected: %s", ', '.join(mapping))

    annotations = []
    gene_ids = set()

    for annotation in gff.parse_gff(input_file):
        annotations.append(annotation)
        gene_ids.add(annotation.gene_id)

    iterator = uniprot_io.uniprot_mappings_to_dict(
        mapping_file,
        gene_ids=set(gene_ids),
        mappings=set(mapping),
        num_lines=None
    )

    if progress:
        iterator = tqdm(iterator)

    file_mappings = dict(iterator)

    count = 0

    for annotation in annotations:
        try:
            mappings = file_mappings[annotation.gene_id]
            for mapping_id, mapping_values in viewitems(mappings):
                if mapping_id == uniprot_io.MAPPINGS['taxonomy']:
                    if (annotation.taxon_id and force_taxon_id) or \
                       (annotation.taxon_id is None):
                        annotation.taxon_id = mapping_values[0]
                        annotation.taxon_db = 'UNIPROT'
                else:
                    annotation.set_mapping(mapping_id, mapping_values)
            count += 1
        except KeyError:
            pass

        annotation.to_file(output_file)

    LOG.info(
        "Number of annotation changed: %d/%d (%.2f%%)",
        count,
        len(annotations),
        count / len(annotations) * 100
    )


def load_htseq_count_files(count_files, samples):
    counts = {}

    for sample, count_file in zip(samples, count_files):
        for uid, count in mgkit.counts.func.load_htseq_counts(count_file):
            if uid not in counts:
                counts[uid] = {}
            counts[uid][sample] = count

    return counts


def load_featurecounts_files(count_files, samples):
    counts = {}
    for sample, count_file in zip(samples, count_files):
        for line in open_file(count_file, 'rb'):
            line = line.decode('ascii')
            if line.startswith('#') or line.lower().startswith('geneid'):
                continue
            line = line.strip().split('\t')
            uid = line[0]
            counts[uid] = dict(zip(samples, line[6:]))
    return counts


@main.command('counts', help="""Adds counts data to the GFF file""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--samples', required=True, multiple=True,
              help="""Sample names, in the same order as the count files""")
@click.option('-c', '--count-files', multiple=True, required=True, help="Count file(s)")
@click.option('-f', '--fpkms', default=False, is_flag=True,
              help="If the counts are FPKMS")
@click.option('-e', '--featureCounts', default=False, is_flag=True,
              help="""If the counts files are from featureCounts""")
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def counts_command(verbose, samples, count_files, fpkms, featurecounts,
                   progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    LOG.info("Samples passed: %s", ','.join(samples))

    if featurecounts:
        counts = load_featurecounts_files(count_files, samples)
    else:
        counts = load_htseq_count_files(count_files, samples)

    key = "fpkms_{}" if fpkms else "counts_{}"

    iterator = gff.parse_gff(input_file)

    if progress:
        iterator = tqdm(iterator)

    for annotation in iterator:
        try:
            ann_counts = counts[annotation.uid]
        except KeyError:
            LOG.warning("No counts found for annotation %s", annotation.uid)
            annotation.to_file(output_file)
            continue

        for sample in samples:
            annotation.attr[key.format(sample)] = ann_counts.get(sample, 0)

        annotation.to_file(output_file)


def parse_hdf5_arg(ctx, param, values):
    if values is None:
        return values
    try:
        file_name, table = values.strip().split(':')
    except ValueError:
        raise click.BadParameter("""The HDF5 file name must be followed by ':'
                                 and the name of the table""")

    return (file_name, table)


@main.command('addtaxa', help='''Adds taxonomy information from a GI-Taxa,
              gene_id/taxon_id table or a dictionary serialised as a
              pickle/msgpack/json file, or a table in a HDF5 file''')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--gene-taxon-table', default=None,
              help="""GIDs taxonomy table (e.g. gi_taxid_nucl.dmp.gz) or a similar file where GENE/TAXON are tab separated and one per line""")
@click.option('-a', '--gene-attr', default='gene_id',
              help="""In which attribute the GENEID in the table is stored (defaults to *gene_id*)""")
@click.option('-f', '--hdf-table', default=None, callback=parse_hdf5_arg,
              help="""HDF5 file and table name to use for taxon_id lookups. The format to pass is the file name, colon and the table file hf5:taxa-table. The index in the table is the accession_id, while the column `taxon_id` stores the taxon_id as int""")
@click.option('-x', '--taxonomy', default=None, help="""Taxonomy file - If given, both *taxon_name* and *lineage* attributes will be set""")
@click.option('-d', '--dictionary', default=None, help="""A serialised dictionary, where the key is the GENEID and the value is TAXONID. It can be in json or msgpack format (can be a compressed file) *Note*: the dictionary values takes precedence over the table files""")
@click.option('-e', '--skip-no-taxon', default=False, is_flag=True,
              help="""If used, annotations with no taxon_id won't be outputted""")
@click.option('-db', '--taxon-db', default='NONE', help="""DB used to add the taxonomic information""")
@click.option('-c', '--cache-table', default=False, is_flag=True,
              help="""If used, annotations are not preloaded, but the taxa table is cached, instead""")
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def addtaxa_command(verbose, gene_taxon_table, hdf_table, gene_attr, taxonomy,
                    dictionary, skip_no_taxon, taxon_db, cache_table,
                    progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    if (gene_taxon_table is not None) and (not cache_table):
        annotations = []
        gene_ids = set()

        # the annotations are preloaded if a table is supplied, such as
        # the one from NCBI, so only the gene_ids necessary are taken from that
        # table to save memory

        for annotation in gff.parse_gff(input_file):
            annotations.append(annotation)
            gene_ids.add(annotation.get_attr(gene_attr, str))
        gene_ids = dict(
            blast.parse_accession_taxa_table(
                gene_taxon_table,
                acc_ids=gene_ids,
                key=0,
                value=1,
                no_zero=True
            )
        )
    elif (gene_taxon_table is not None) and cache_table:
        annotations = gff.parse_gff(input_file)
        gene_ids = cache_dict_file(
            blast.parse_accession_taxa_table(
                gene_taxon_table,
                key=0,
                value=1,
                no_zero=True
            )
        )
    elif hdf_table is not None:
        try:
            gene_ids = HDFDict(hdf_table[0], hdf_table[1])
        except ValueError as e:
            utils.exit_script(
                str(e),
                3
            )
        annotations = gff.parse_gff(input_file)
    else:
        # in case a dictionary is supplied, it's expected to fit in memory,
        # meaning that the GFF doesn't have to be preloaded
        annotations = gff.parse_gff(input_file)
        gene_ids = {}

    if (dictionary is not None) and (hdf_table is None):
        dict_file = open_file(dictionary, 'r')
        if '.json' in dictionary:
            gene_ids.update(json.load(dict_file))
        elif '.msgpack':
            gene_ids.update(msgpack.load(dict_file))
        else:
            gene_ids.update(pickle.load(dict_file))

        # Ensures all taxon_ids are *int*
        gene_ids = dict(
            (key, int(value))
            for key, value in viewitems(gene_ids)
        )

    if taxonomy is not None:
        taxonomy = taxon.Taxonomy(taxonomy)

    if progress:
        annotations = tqdm(annotations)

    count = 0
    for lineno, annotation in enumerate(annotations):
        gene_id = annotation.get_attr(gene_attr, str)
        try:
            taxon_id = gene_ids[gene_id]
            annotation.taxon_id = taxon_id
            annotation.set_attr('taxon_db', taxon_db)
            count += 1
        except KeyError:
            LOG.warning(
                "No Taxon ID for GENE - %s",
                gene_id
            )
            taxon_id = None
            if skip_no_taxon:
                continue

        if (taxonomy is not None) and (taxon_id is not None):
            try:
                annotation.attr['taxon_name'] = taxonomy[taxon_id].s_name
                annotation.attr['lineage'] = ','.join(
                    taxon_name
                    for taxon_name in taxon.get_lineage(taxonomy, taxon_id, names=True, only_ranked=True)
                    if taxon_name
                )
            except KeyError:
                LOG.warning("Taxon ID %d not found in the Taxonomy", taxon_id)

        annotation.to_file(output_file)

    LOG.info("Percentage of annotations changed %.2f%%", count / (count + 1) * 100)


@main.command('pfam', help="""Adds information from Pfam""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-i', '--id-attr', default='gene_id', help="""In which attribute the Pfam ID/ACCESSION is stored (defaults to *gene_id*)""")
@click.option('-a', '--use-accession', default=False, is_flag=True, help="""If used, the attribute value is the Pfam ACCESSION (e.g. PF06894), not ID (e.g. Phage_TAC_2)""")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def pfam_command(verbose, id_attr, use_accession, input_file, output_file):
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    LOG.info("Downloading Pfam Information")

    pfam_families = pfam.get_pfam_families(
        'acc' if use_accession else 'id'
    )

    for annotation in gff.parse_gff(input_file):
        pfam_id = annotation.attr[id_attr]
        try:
            annotation.attr['pfam_description'] = pfam_families[pfam_id][1]
        except KeyError:
            LOG.warning("No description found for family %s", pfam_id)
        annotation.to_file(output_file)


def update_annotation_coverage(annotations, depth, sample):
    for annotation in annotations:
        cov = depth.region_coverage(
            annotation.seq_id,
            annotation.start,
            annotation.end
        )
        if sample is None:
            annotation.set_attr('cov', cov)
        else:
            annotation.set_attr(sample, cov)



@main.command('cov_samtools', help="""Adds information from samtools_depth""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('-m', '--average', is_flag=True, help='if one or more samples are provided, the average coverage is calculated')
@click.option('-s', '--samples', default=None, multiple=True,
              help='''Sample name, will add a `sample_cov` in the GFF file. If not passed, the attribute will be `cov`''')
@click.option('-d', '--depths', required=True, multiple=True,
              help='`samtools depth -aa` file')
@click.option('-n', '--num-seqs', default=0,  type=click.INT,
              show_default=True, help='''Number of sequences to update the log. If 0, no message is logged''')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def samtools_depth_command(verbose, average, samples, depths, num_seqs, progress,
                           input_file, output_file):
    logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if not samples:
        samples = (None,)
    else:
        samples = ['{}_cov'.format(sample) for sample in samples]

    max_size_dict = {}
    annotations = gff.group_annotations(gff.parse_gff(input_file), key_func=lambda x: x.seq_id)

    if progress:
        it = tqdm(annotations.items(), desc='Max lengths')
    else:
        it = annotations.items()
    for seq_id, seq_annotations in it:
        max_size_dict[seq_id] = max(annotation.end for annotation in seq_annotations)

    depths = [
        align.SamtoolsDepth(
            file_name,
            num_seqs=None if num_seqs == 0 else num_seqs,
            max_size_dict=max_size_dict,
        )
        for file_name in depths
    ]
    if len(samples) != len(depths):
        utils.exit_script('Number of samples different from number of files', 2)

    if progress and (len(samples) > 1):
        depth_it = tqdm(zip(samples, depths), desc='Samples', total=len(samples))
    else:
        depth_it = zip(samples, depths)

    for sample, depth in depth_it:
        seq_found = set()

        if progress:
            seq_bar = tqdm(desc='Sequences', total=len(annotations))
        while True:
            seq_id = depth.advance_file()
            # sequence is in the annotations
            if seq_id in annotations:
                update_annotation_coverage(annotations[seq_id], depth, sample)
                seq_found.add(seq_id)
                depth.drop_sequence(seq_id)
                if progress:
                    seq_bar.update(n=1)
            # reached the end of the file
            # updated the remaining sequences
            elif seq_id is None:
                for not_found in set(annotations) - seq_found:
                    update_annotation_coverage(annotations[not_found], depth,
                                                sample)
                    if progress:
                        seq_bar.update(n=1)
                LOG.warning('Cannot find %d of %d sequences in sample %s depth file',
                    len(depth.not_found), len(annotations), sample)
                break
            # found all sequences, no need to continue scanning the file
            elif len(seq_found) == len(annotations):
                LOG.info("Found all sequences in GFF file (sample %s)", sample)
                break
            # the sequence was not in the annotations, skip
            # continue scanning the depth file and drop the seq_id
            # (redundat?)
            elif seq_id not in annotations:
                depth.drop_sequence(seq_id)

            # some debug code, to check the density of the SparseArray(s)
            # if (num_seqs > 0) and (len(depth.density) % num_seqs == 0):
            #     LOG.debug(sum(depth.density) / len(depth.density))

        # closes the progress bar
        if progress:
            seq_bar.close()

    # If no sample was provided it's assumed that the average was calculated
    if average and (samples[0] is not None):
        max_coverage = 0

        LOG.info('Calculating average coverage')
        for seq_annotations in annotations.values():
            for annotation in seq_annotations:
                d = annotation.sample_coverage
                annotation.set_attr('cov', sum(d.values()) / len(d))
                max_coverage = max(max_coverage, max(d.values()))
        LOG.info('Maximum Coverage found: %d', max_coverage)

    LOG.info('Writing annotations to file (%r)', output_file)
    for annotation in itertools.chain(*annotations.values()):
        annotation.to_file(output_file)
