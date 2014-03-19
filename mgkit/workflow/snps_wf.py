"""
Workflow associated with SNPs
"""

import cPickle
import logging
import numpy
# from ..unused.snps_old import combine_genes_by_feature, order_ratios
from ..mappings.utils import count_genes_in_mapping
# from ..plots import boxplot_snp, barchart_categories
from ..utils.dictionary import split_dictionary_by_value, reverse_mapping

LOG = logging.getLogger(__name__)


def load_sample_data(sample_data):
    LOG.info("Loading sample data from file %s", sample_data.name)
    data = cPickle.load(sample_data)
    return data


def write_single_values(values, sign_values, id_name, out_file, sep='\t'):
    LOG.info("Writing results to file")
    out_file = open(out_file, 'w')

    out_file.write(
        sep.join(
            (
                id_name,
                "significant",
                "pvalue",
                "mean pN/pS",
                "num replicates",
                "all pN/pS"
            )
        ) + '\n'
    )
    for value_id, ratios in values.iteritems():
        issign = 'Y' if value_id in sign_values else 'N'
        try:
            pvalue = sign_values[value_id]
        except KeyError:
            pvalue = ''
        out_file.write(
            "{value_id}\t{sign}\t{pvalue}\t{meanp}\t{nratios}\t{ratios}\n".
            format(
                value_id=value_id,
                sign=issign,
                pvalue=pvalue,
                meanp=numpy.mean(ratios),
                nratios=len(ratios),
                ratios=sep.join(str(ratio) for ratio in ratios)
            )
        )


def get_sample_ratios(sample_data, min_cov):
    s_values = combine_genes_by_feature(
        sample_data,
        'sample',
        min_cov
    )
    return [sample.calc_ratio() for sample in s_values.itervalues()]


def write_gt_values(gt_values, sign_genes, out_file, dnds=None, sep='\t'):
    LOG.info("Writing Gene-Taxon values to file")
    # if dnds is not None:
    #     LOG.info("dN/dS data available")
    out_file = open(out_file, 'w')

    out_file.write(
        sep.join(
            (
                "gene",
                "taxon",
                "significant",
                "higher",
                "mean dN/dS",
                "mean pN/pS",
                "num replicates",
                "all pN/pS"
            )
        ) + '\n'
    )

    for gene_id, taxa_dict in gt_values.iteritems():
        issign = 'Y' if gene_id in sign_genes else 'N'
        for taxon_name, ratios in taxa_dict.iteritems():
            if (issign == 'Y') and (taxon_name in sign_genes[gene_id]):
                ishigher = 'Y'
            else:
                ishigher = 'N'
            if dnds is not None:
                meand = dnds.get((gene_id, taxon_name), '')
            else:
                meand = ''
            out_file.write(
                (
                    "{gene_id}\t{taxon_name}\t{sign}\t{higher}\t" +
                    "{meand}\t{meanp}\t{nratios}\t{ratios}\n"
                ).format(
                    gene_id=gene_id,
                    taxon_name=taxon_name,
                    sign=issign,
                    higher=ishigher,
                    meand=meand,
                    meanp=numpy.mean(ratios),
                    nratios=len(ratios),
                    ratios=sep.join(str(ratio) for ratio in ratios)
                )
            )


def plot_results_boxplot(id_name, s_values, values, sign_values, fontsize,
                         aspect, key_filter, base_name, ext='.pdf', colors=None):
    """
    .. todo::

        move split dictionary in another function?
    """
    LOG.info("Saving all %s plot", id_name)
    plot_order = order_ratios(values)
    plot_file = base_name + 'all' + ext
    boxplot_snp(values, plot_order, labelfont=fontsize, file_name=plot_file,
                title="pN/pS among samples for all " + id_name, ylabel='pN/pS',
                xlabel=id_name.capitalize(), fig_aspect=aspect,
                taxon_colors=colors)
    sign_genes_low, sign_genes_high = split_dictionary_by_value(
        values,
        numpy.median(s_values),
        key_filter=key_filter
    )
    plot_file = base_name + 'sign' + ext
    LOG.info("Saving significant %s plots", id_name)
    plot_order = order_ratios(values, key_filter=key_filter)
    boxplot_snp(values, plot_order, labelfont=fontsize, file_name=plot_file,
                title="pN/pS among samples for significant " + id_name,
                ylabel='pN/pS', xlabel=id_name.capitalize(), fig_aspect=aspect,
                taxon_colors=colors)
    plot_file = base_name + 'sign-low' + ext
    plot_order = order_ratios(sign_genes_low)
    boxplot_snp(values, plot_order, labelfont=fontsize, file_name=plot_file,
                title="pN/pS among samples for significant {0}\n".
                format(id_name) + "(lower than samples median)",
                ylabel='pN/pS', xlabel=id_name.capitalize(), fig_aspect=aspect,
                taxon_colors=colors)
    plot_file = base_name + 'sign-high' + ext
    plot_order = order_ratios(sign_genes_high)
    boxplot_snp(values, plot_order, labelfont=fontsize, file_name=plot_file,
                title="pN/pS among samples for significant {0}\n".
                format(id_name) + "(higher than samples median)",
                ylabel='pN/pS', xlabel=id_name.capitalize(), fig_aspect=aspect,
                taxon_colors=colors)


def plot_gene_taxon_barcharts(sign_high, mappings, names, labels):
    #plots barcharts for gene-taxon dictionary
    # mappings = (egmap, gos_gen_map, gos_met_map, cz_map)
    # names = (eggnog.EGGNOG_CAT, gos_gen_names, gos_met_names, None)
    # labels = ('eggNOG', 'GOSlim Generic', 'GOSlim Metagenomic', 'CaZy')

    taxa = set(reverse_mapping(sign_high).keys())

    for taxon_name in taxa:
        low_genes = sign_low_rev.get(taxon_name, [])
        tot_genes = reverse_mapping(gt_vals).get(taxon_name, [])
        for mapping, name, label in zip(mappings, names, labels):
            print taxon_name, label

            try:
                df = count_genes_in_mapping(
                    (low_genes, tot_genes),# sign_high),
                    ('Higher than median', 'Total'),# 'Total significant'),
                    mapping, normalise=True
                )
            except ZeroDivisionError:
                print "\tSkip", taxon_name, label
                continue
            if len(df.columns) == 0:
                print "\tSkip", taxon_name, label
                continue
            barchart_categories(df, colors=('g','y', 'm'), xlabel_dict=name, fig_dpi=300,
                           title="{} - {}".format(taxon_name, label), fig_size=(0.5 * len(df.columns), 15.),
                           file_name="cat_counts/{}_{}.pdf".format(taxon_name, label))


def plot_results_categories(values, sign_low, sign_high, fontsize,
                            mappings, names, labels,
                            aspect, base_name, ext='.pdf', colors=None):
    for label, mapping, name in zip(labels, mappings, names):
        table = count_genes_in_mapping(
            (sign_low, sign_low, values)
            ('Lower than median', 'Higher than median', 'Total'),
            mapping, normalise=True
        )
        barchart_categories(
            table, colors=colors, xlabel_dict=name, fig_aspect=aspect,
            title="Gene - " + label, file_name=base_name + label + ext)
