"""
Add more information to GFF annotations: gene mappings, coverage, taxonomy,
etc..

Uniprot Command
***************

If the *gene_id* of an annotation is a Uniprot ID, the script queries Uniprot
for the requested information. At the moment the information that can be added
is the taxon_id, taxon_name, lineage and mapping to EC, KO, eggNOG IDs.

It's also possible to add mappings to other databases using the *-m* option with
the correct identifier for the mapping, which can be found at `this page
<http://www.uniprot.org/faq/28>`_; for example if it's we want to add the
mappings of uniprot IDs to *BioCyc*, in the *abbreviation* column of the
mappings we find that it's identifier is *REACTOME_ID*, so we pass *-m REACTOME*
to the script (leaving *_ID* out). Mapped IDs are separated by commas.

The taxonomy IDs are not overwritten if they are found in the annotations, the
*-f* is provided to force the overwriting of those values.

See also :ref:`gff-specs` for more informations about the GFF specifications
used.

.. note::

    As the script needs to query Uniprot a lot, it is reccommended to split
    the GFF in several files, so an error in the connection doesn't waste time.

Taxonomy Command
****************

To refine the taxonomic assignments of predicted genes annotations, the
annotation sequences may be searched against a database like the NCBI *nt*.

This commands takes as input a GFF file, one or more blast output files and a
file with all mappings from GIDs to taxonomy IDs. More information on how to
get the file can be read in the documentation of the function
:func:`mgkit.io.blast.parse_gi_taxa_table`.

The fasta sequences used with BLAST must have as name the uid of the annotations
they refer to, and one way to obtain these sequences is to use the function
:func:`mgkit.io.gff.extract_nuc_seqs` and save them to a fasta file. Another
options is to use the `sequence` command of the `get-gff-info` script
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
`sample,samplealg.bam`. More samples can be supplied adding more `-a` arguments.

.. hint::

    As an example, to add coverage for `sample1`, `sample2` the command line
    is::

        add-gff-info coverage -a sample1,sample1.bam -a sample2,sample2.bam \\
        inputgff outputgff

A total coverage for the annotation is also calculated and stored in the
`cov` attribute, while each sample coverage is stored into `sample_cov` as per
:ref:`gff-specs`.

Uniprot Offline Mappings
************************

Similar to the *uniprot* command, it uses the `idmapping <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_
file provided by Uniprot, which speeds up the process of adding mappings and
taxonomy IDs from Uniprot gene IDs. It's not possible tough to add *EC* mappings
with this command, as those are not included in the file.

Expected Aminoacidic Changes
****************************

Some scripts, like :ref:`snp-parser`, require information about the expected
number of synonymous and non-synonymous changes of an annotation. This can be
done using :meth:`mgkit.io.gff.Annotation.add_exp_syn_count` by the user of the
command `exp_syn` of this script. The attributes added to each annotation are
explained in the :ref:`gff-specs`

Changes
*******

.. versionadded:: 0.1.12

.. versionchanged:: 0.1.13

* added *--force-taxon-id* option to the *uniprot* command
* added *coverage* command
* added *taxonomy* command
* added *unipfile* command

.. versionchanged:: 0.1.15
    *taxonomy* command *-b* option changed

.. versionchanged:: 0.1.16
    added *exp_syn* command

"""
from __future__ import division
import sys
import argparse
import logging
import itertools
import functools
import pysam
import mgkit
from . import utils
from .. import align
from .. import logger
from .. import taxon
from ..io import gff, blast, fasta
from ..io import uniprot as uniprot_io
from ..net import uniprot as uniprot_net

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
    group = parser.add_argument_group('Require Internet connection')
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
        '-m',
        '--mapping',
        action='append',
        type=str,
        help='Add any DB mappings to annotations'
    )

    parser.set_defaults(func=uniprot_command)


def add_uniprot_info(annotations, options):
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
        columns.append('lineage()')
    if options.mapping is not None:
        for db in options.mapping:
            columns.append('database({0})'.format(db))

    if not columns:
        return

    LOG.info("Retrieving gene information from Uniprot")

    data = uniprot_net.get_gene_info(
        [x.gene_id for x in annotations],
        columns=columns,
        contact=options.email
    )

    for annotation in annotations:
        try:
            gene_info = data[annotation.gene_id]
        except KeyError:
            #no data was found
            continue

        for column, values in gene_info.iteritems():
            #nothing found
            if not values:
                continue

            if column == 'organism-id':
                if (annotation.taxon_id and options.force_taxon_id) or \
                   (annotation.taxon_id is None):
                    annotation.attr['taxon_id'] = int(values)
                    annotation.attr['taxon_db'] = 'UNIPROT'
                    #test with a try/expect maybe
                    if 'organism' in columns:
                        annotation.attr['taxon_name'] = gene_info['organism']
                    if column.startswith('lineage'):
                        annotation.attr['lineage'] = gene_info['lineage()']
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


def uniprot_command(options):

    if options.buffer < 1:
        options.buffer = 1

    ann_buffer = []

    for annotation in gff.parse_gff(options.input_file, gff_type=gff.from_gff):

        ann_buffer.append(annotation)

        if len(ann_buffer) == options.buffer:

            add_uniprot_info(ann_buffer, options)

            for annotation in ann_buffer:
                annotation.to_file(options.output_file)

            ann_buffer = []
    else:
        add_uniprot_info(ann_buffer, options)

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
    #the minimum score required to be part of the taxon IDs selected is to
    #be at most 10 bits (by default) from the maximum bitscore found in the hits
    min_score = max(hits, key=lambda x: x[-1])[-1] - score

    #filter over only cellular organisms
    cellular_organism = 131567

    taxon_ids = set()

    for gid, identity, bitscore in hits:
        #No taxon_id found for the gid
        taxon_id = gid_taxon_map.get(gid, None)

        #skip already included taxon_ids
        if taxon_id in taxon_ids:
            continue

        #if the requirements are met
        if (taxon_id is not None) and (bitscore >= min_score):
            try:
                #test to make sure it is a cellular organism
                if taxon.is_ancestor(taxonomy, taxon_id, cellular_organism):
                    taxon_ids.add(taxon_id)
                #not part of the loop
                else:
                    LOG.debug("%d not part of cellular organisms", taxon_id)
            #The taxon id is not in the taxonomy used. It's skipped
            except KeyError:
                LOG.warning("%d is not part of the taxonomy", taxon_id)
                continue

    #no taxon_id passes the filters
    if len(taxon_ids) == 0:
        return None

    #log message fired up only if the package is set on DEBUG
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
    parser.set_defaults(func=exp_syn_command)


def exp_syn_command(options):
    """
    .. versionadded:: 0.1.16
    """

    seqs = dict(fasta.load_fasta(options.reference))

    for annotation in gff.parse_gff(options.input_file):
        annotation.add_exp_syn_count(seqs[annotation.seq_id])

        annotation.to_file(options.output_file)


def uniprot_offline_command(options):

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
    utils.add_basic_options(parser_u)

    parser_t = subparsers.add_parser(
        'taxonomy',
        help='''Adds taxonomic information from annotation sequences blasted
                against a NCBI db'''
    )

    set_blast_taxonomy_parser(parser_t)
    set_common_options(parser_t)
    utils.add_basic_options(parser_t)

    parser_c = subparsers.add_parser(
        'coverage',
        help='Adds coverage information from BAM Alignment files'
    )

    set_coverage_parser(parser_c)
    set_common_options(parser_c)
    utils.add_basic_options(parser_c)

    parser_e = subparsers.add_parser(
        'exp_syn',
        help='Adds expected synonymous and non-synonymous changes information'
    )

    set_exp_syn_parser(parser_e)
    set_common_options(parser_e)
    utils.add_basic_options(parser_e)

    parser_f = subparsers.add_parser(
        'unipfile',
        help='Adds mappings and taxonomy from Uniprot mapping file'
    )

    set_uniprot_offline_parser(parser_f)
    set_common_options(parser_f)
    utils.add_basic_options(parser_f)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)
