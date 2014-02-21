#!/usr/bin/env python
"""

passare in summarise_gff con relative tabelle?
    coverage_of_taxa-min10
    number_of_taxa_per_gene
    number_of_genes_per_taxon

OK:
    fare qui
    taxa-all-sorted_by_median-min10
    boxplot-all_genes-sorted_by_median-min10

    wilcoxon test:
        geni
        taxa
--------------------------------------------

wilcoxon pairwise test:
    boxplot per taxon con tutti i geni significativi per quel taxon (se e' uno
    solo il gene evitare tight layout)
    barchart delle categorie per taxon

clustering:
    usare PCA per decidere le componenti in automatico? riservare per
    successivo sviluppo o lasciare primo quartile e media come misure?
    - non fare per il momento. Se fatto, meglio usare le tre misure (primo
    quartile, media e secondo quartile) ed un plot 3D

    i dati dnds derivano da tutti i geni hanno almeno 10 repliche che superano
    il test wilcoxon pairwise (fatto)
    poter scegliere se usare solo quelli significativi dal test wilcoxon o meno
    workflow dN/dS: saltare profile con solo una sequenza - fare successivamente

    costruire script per generare in automatico i valori di dN/dS con crann
        - generare il plot
        - generare la lista dei geni-taxa e loro appartenenza al cluster in un
          foglio excel

term enrichment:
    integrare goatools?
    ouput dei dati comunque nel formato richiesto
    controllare invece possibili altri tool usati da altre persone

"""

import cPickle
import argparse
import logging
import os
import numpy
import mgkit
from mgkit import taxon
from mgkit import logger
from mgkit.taxon import BLACK_LIST
from mgkit.workflow.snps_wf import load_sample_data, write_single_values
from mgkit.workflow.snps_wf import write_gt_values, plot_results_categories
from mgkit.workflow.snps_wf import get_sample_ratios, plot_results_boxplot
from mgkit.snps import combine_genes_by_feature
from mgkit.snps import convert_dict_for_pairwise_wilcoxon
from mgkit.snps import wilcoxon_test_ratio_pairwise
from mgkit.snps import convert_dict_for_wilcoxon, wilcoxon_test_ratio
from mgkit.plots import get_taxon_colors
from mgkit.utils.dictionary import split_dictionary_by_value

#if r_utils modules cannot be imported gives a warning
try:
    from mgkit.utils.r_func import correct_pvalues
except ImportError:
    correct_pvalues = None


LOG = logging.getLogger(__name__)


def make_dir(dir_name):
    try:
        LOG.info("Writing ouput files to %s", dir_name)
        os.mkdir(dir_name)
    except OSError:
        LOG.warning("Output directory %s already present", dir_name)


def get_taxa_colors(taxa, taxon_data):
    root_map = taxon.load_taxonomy_map(taxon_data)
    color_map = get_taxon_colors(taxa, root_map)
    return color_map


def load_mapping_data(map_files):
    pass


def gene_taxon_clustering(gt_values):
    pass


def taxon_analysis(options):
    LOG.info("Starting Taxon analysis")
    sample_data = load_sample_data(options.sample_data)
    make_dir(options.out_dir)
    t_values = combine_genes_by_feature(
        sample_data,
        'taxon',
        options.min_cov
    )
    s_values = get_sample_ratios(sample_data, options.min_cov)
    t_values = convert_dict_for_wilcoxon(t_values)
    #distribuzione campio
    sign_taxa = wilcoxon_test_ratio(
        s_values,
        t_values,
        corr_func=correct_pvalues if options.no_corr else None,
        corr_method=options.corr_method,
        min_num=options.min_num,
        threshold=options.threshold
    )
    LOG.info("Number of significant taxa: %d", len(sign_taxa))
    write_single_values(
        t_values,
        sign_taxa,
        'taxon',
        os.path.join(
            options.out_dir,
            'taxon_values.tab'
        )
    )
    color_map = get_taxa_colors(t_values, options.taxon_data)
    plot_base = os.path.join(options.out_dir, 'taxa_plot-')
    plot_results_boxplot('taxa', s_values, t_values, sign_taxa,
                         options.label_size, options.fig_ratio,
                         sign_taxa.index, plot_base, colors=color_map)


def gene_analysis(options):
    LOG.info("Starting Gene analysis")
    sample_data = load_sample_data(options.sample_data)
    make_dir(options.out_dir)
    g_values = combine_genes_by_feature(
        sample_data,
        'gene',
        options.min_cov
    )
    g_values = convert_dict_for_wilcoxon(g_values)
    s_values = get_sample_ratios(sample_data, options.min_cov)
    #distribuzione campio
    sign_genes = wilcoxon_test_ratio(
        s_values,
        g_values,
        corr_func=correct_pvalues if options.no_corr else None,
        corr_method=options.corr_method,
        min_num=options.min_num,
        threshold=options.threshold
    )
    LOG.info("Number of significant genes: %d", len(sign_genes))
    write_single_values(
        g_values,
        sign_genes,
        'gene',
        os.path.join(
            options.out_dir,
            'gene_values.tab'
        )
    )
    plot_base = os.path.join(options.out_dir, 'genes_plot-')
    plot_results_boxplot('genes', s_values, g_values, sign_genes,
                         options.label_size, options.fig_ratio,
                         sign_genes.index, plot_base)
    sign_genes_low, sign_genes_high = split_dictionary_by_value(
        g_values,
        numpy.median(s_values),
        key_filter=sign_genes.index
    )
    plot_base = os.path.join(options.out_dir, 'genes_categories-')
    plot_results_categories(g_values,
                            sign_genes_low,
                            sign_genes_high,

                            options.label_size,
                            aspect=options.fig_ratio,
                            colors=('g', 'r', 'b'),
                            base_name=plot_base
                            )


def gene_taxon_analysis(options):
    LOG.info("Starting Gene-Taxon analysis")
    sample_data = load_sample_data(options.sample_data)
    make_dir(options.out_dir)
    gt_values = combine_genes_by_feature(
        sample_data,
        'gene-taxon',
        options.min_cov
    )
    gt_values = convert_dict_for_pairwise_wilcoxon(gt_values)
    sign_genes = wilcoxon_test_ratio_pairwise(
        gt_values,
        corr_func=correct_pvalues if options.no_corr else None,
        corr_method=options.corr_method,
        min_num=options.min_num,
        threshold=options.threshold
    )
    LOG.info("Number of significant genes: %d", len(sign_genes))
    dnds = None
    if options.dnds_file is not None:
        LOG.info("Loading dN/dS data from file %s", options.dnds_file.name)
        dnds = cPickle.load(options.dnds_file)
        # gene_taxon_clustering(gt_values)

    write_gt_values(
        gt_values, sign_genes, os.path.join(
            options.out_dir,
            'gene-taxon_values.tab'
        ),
        dnds=dnds
    )


def add_stats_options(parser):
    parser.add_argument(
        '-c',
        '--no-corr',
        default=True,
        action='store_false',
        help='Disables p-value correction using R'
    )
    parser.add_argument(
        '-t',
        '--corr-method',
        default='bonferroni',
        action='store',
        help='Method used for p-value correction: bonferroni, BH, etc.'
    )
    parser.add_argument(
        '--min-num',
        default=10,
        action='store',
        type=int,
        help='minimum number of replicates accepted for a gene/taxon'
    )
    parser.add_argument(
        '--threshold',
        default=0.05,
        action='store',
        type=int,
        help='p-value threshold for acceptance'
    )


def add_graph_options(parser):
    """
    Sets command line arguments parser
    """
    parser.add_argument(
        '--fig-ratio',
        default=0.3,
        action='store',
        type=float,
        help='Figure aspect ratio'
    )
    parser.add_argument(
        '--label-size',
        default=6,
        action='store',
        type=int,
        help='X axis ticks font size. 0 means no labels'
    )


def add_general_options(parser):
    """
    Sets command line arguments parser
    """
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_const',
        const=logging.DEBUG,
        default=logging.INFO,
        help='more verbose'
    )
    parser.add_argument(
        '-m',
        '--min-cov',
        default=4,
        action='store',
        type=int,
        help='Minimum per sample coverage accepted'
    )
    parser.add_argument(
        '-b',
        '--black-list',
        nargs='*',
        default=BLACK_LIST,
        # action='store',
        help='Taxa to be excluded from analysis'
    )


def add_taxon_command(subparser):
    parser = subparser.add_parser(
        'taxon',
        help='Taxon analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--taxon-data',
        default='data/taxonomy.pickle',
        action='store',
        help='Taxonomy data, used to color plots'
    )
    parser.set_defaults(func=taxon_analysis)
    add_general_options(parser)
    add_stats_options(parser)
    add_command_commons(parser)


def add_gene_command(subparser):
    parser = subparser.add_parser(
        'gene',
        help='Gene analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--eggnog-data',
        default=None,
        action='store',
        help='eggNOG data, used for mappings'
    )
    parser.add_argument(
        '--cazy-data',
        default=None,
        action='store',
        help='CaZy data, used for mappings'
    )
    parser.add_argument(
        '-g',
        '--go-data',
        default=None,
        action='store',
        help='GO data, used for mappings'
    )
    parser.set_defaults(func=gene_analysis)
    add_general_options(parser)
    add_stats_options(parser)
    add_command_commons(parser)


def add_gene_taxon_command(subparser):
    parser = subparser.add_parser(
        'gene-taxon',
        help='Gene-Taxon analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.set_defaults(func=gene_taxon_analysis)
    parser.add_argument(
        '-d',
        '--dnds-file',
        type=argparse.FileType('r'),
        action='store',
        default=None,
        help='File containing data for dN/dS values'
    )
    add_general_options(parser)
    add_stats_options(parser)
    add_command_commons(parser)


def add_command_commons(parser):
    parser.add_argument(
        'sample_data',
        type=argparse.FileType('r'),
        action='store',
        help='Saved data from VCF and SNPDat parsing'
    )
    parser.add_argument(
        'out_dir',
        action='store',
        help='Directory to which results are saved'
    )
    add_graph_options(parser)


def main():
    parser = argparse.ArgumentParser(
        description='SNPs analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {0}'.format(mgkit.__VERSION__)
    )

    subparser = parser.add_subparsers(help='Possible analysis')
    add_gene_taxon_command(subparser)
    add_taxon_command(subparser)
    add_gene_command(subparser)

    options = parser.parse_args()

    logger.config_log(options.verbose)

    if correct_pvalues is None:
        LOG.warning(
            "Error in importing r_utils, p-value correction won't be performed"
        )

    options.func(options)

if __name__ == '__main__':
    main()
