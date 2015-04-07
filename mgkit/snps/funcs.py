"""
Functions used in SNPs manipulation
"""
import logging
import numpy
import itertools
import functools
import pandas
import scipy.stats
import csv
import copy
from .classes import GeneSyn
from .filter import pipe_filters

LOG = logging.getLogger(__name__)


def build_rank_matrix(dataframe, taxonomy=None, taxon_rank=None):
    """
    Make a rank matrix from a :class:`pandas.Series` with the pN/pS values of a
    dataset.

    :param dataframe: :class:`pandas.Series` instance with a MultiIndex
        (gene-taxon)
    :param taxonomy: :class:`taxon.UniprotTaxonomy` instance with the full
        taxonomy
    :param str taxon_rank: taxon rank to limit the specifity of the taxa
        included

    :return: :class:`pandas.DataFrame` instance
    """
    #aggiungere controllo che sia un vettore, altrimenti convertire
    if taxon_rank is None:
        taxon_index = sorted(set(dataframe.index.get_level_values('taxon')))
    else:
        taxon_index = sorted(
            set(
                taxonomy.get_ranked_taxon(taxon_id, taxon_rank).taxon_id
                for taxon_id in set(dataframe.index.get_level_values('taxon'))
                if taxonomy.get_ranked_taxon(
                    taxon_id, taxon_rank
                ).rank == taxon_rank
            )
        )

    gene_index = sorted(set(dataframe.index.get_level_values('gene')))

    # print taxon_index
    rank_matrix = pandas.DataFrame(index=gene_index, columns=taxon_index)

    for gene_id in rank_matrix.index:
        gene_array = dataframe.loc[gene_id]
        # print gene_id, type(gene_array)
        iterator = zip(gene_array.index, scipy.stats.rankdata(gene_array))
        for taxon_id, rank in iterator:
            if taxon_rank is not None:
                taxon_id = taxonomy.get_ranked_taxon(
                    taxon_id, taxon_rank
                ).taxon_id
                if taxonomy[taxon_id].rank != taxon_rank:
                    continue

            old_value = rank_matrix.get_value(gene_id, taxon_id)
            if numpy.isnan(old_value) or old_value > rank:
                rank_matrix.set_value(gene_id, taxon_id, rank)

    return rank_matrix


def group_rank_matrix(dataframe, gene_map):
    """
    Group a rank matrix using a mapping, in the form map_id->ko_ids.

    :param dataframe: instance of a rank matrix from :func:`build_rank_matrix`
    :param dict gene_map: dictionary with the mapping

    :return: :class:`pandas.DataFrame` instance
    """
    rank_matrix = pandas.DataFrame(
        index=gene_map.keys(),
        columns=dataframe.columns
    )

    for mapping_id, gene_ids in gene_map.iteritems():
        mapped_matrix = dataframe.loc[gene_ids]
        #we only use the minimum rank among the genes with a set function
        #min() will return only those
        for taxon_id, rank in mapped_matrix.mean().dropna().iteritems():
            rank_matrix.set_value(mapping_id, taxon_id, rank)

    return rank_matrix


def write_sign_genes_table(out_file, dataframe, sign_genes, taxonomy,
                           gene_names=None):
    """
    Write a table with the list of significant genes found in a dataframe, the
    significant gene list is the result of
    :func:`wilcoxon_pairwise_test_dataframe`.

    :out_file: the file name or file object to write the file
    :dataframe: the dataframe which was tested for significant genes
    :sign_genes: gene list that are significant
    :taxonomy: taxonomy object
    :gene_names: dictionary with the name of the the genes. Optional
    """
    if isinstance(out_file, str):
        out_file = open(out_file, 'w')

    out_file = csv.writer(out_file)

    taxon_ids = sorted(set(dataframe.index.get_level_values('taxon')),
                       key=lambda x: taxonomy[x].s_name)

    genes = set(dataframe.index.get_level_values('gene'))

    header = ['gene_id', 'gene_name', 'significant', 'number of taxa']

    for taxon_id in taxon_ids:
        header.append(taxonomy[taxon_id].s_name + '_mean')
        header.append(taxonomy[taxon_id].s_name + '_stdev')

    out_file.writerow(header)

    for gene in genes:
        values = [''] * (len(taxon_ids) * 2)

        dataframe_gene = dataframe.loc[gene]

        for taxon_id, series in dataframe_gene.iterrows():
            value_idx = taxon_ids.index(taxon_id) * 2

            values[value_idx] = series.mean()
            values[value_idx + 1] = series.std()

        out_file.writerow(
            [
                gene,
                '' if gene_names is None else gene_names[gene],
                'Y' if gene in sign_genes else 'N',
                len(dataframe_gene)
            ] + values
        )


def order_ratios(ratios, aggr_func=numpy.median, reverse=False,
                 key_filter=None):
    """
    Given a dictionary of id->iterable where iterable contains the values of
    interest, the function uses aggr_func to sort (ascending by default) it and
    return a list with the key in the sorted order.

    :param dict ratios: dictionary instance id->iterable
    :param function aggr_func: any function returning a value that can be used
        as a key in sorting
    :param bool reverse: the default is ascending sorting (False), set to True
        to reverse key_filter: list of keys to use for ordering, if None, every
        key is used

    :return: iterable with the sort order
    """

    if key_filter is None:
        key_filter = ratios.keys()

    order = [
        (
            aggr_func(
                [ratio for ratio in ratios[key] if not numpy.isnan(ratio)]
            ) if [ratio for ratio in ratios[key] if not numpy.isnan(ratio)]
            else 0,
            key
        )
        for key in key_filter
    ]

    order.sort(key=lambda x: x[0], reverse=reverse)
    # print order

    # for x, y in order: print y, x, ratios[y]

    return [x[1] for x in order]


def combine_sample_snps(snps_data, min_num, filters, index_type=None,
                        gene_func=None, taxon_func=None):
    """
    Combine a dictionary sample->gene_index->GeneSyn into a
    :class:`pandas.DataFrame`. The dictionary is first filtered with the
    functions in `filters`, mapped to different taxa and genes using
    `taxon_func` and `gene_func` respectively. The returned DataFrame is also
    filtered for each row having at least a `min_num` of not NaN values.

    .. todo::

        detail usage and examples.

    Arguments:
        snps_data (dict): dictionary with the `GeneSNP` instances
        min_num (int): the minimum number of not NaN values necessary in a row
            to be returned
        filters (iterable): iterable containing filter functions, a list can be
            found in :mod:`mgkit.snps.filter`
        index_type (str, None): if `None`, each row index for the DataFrame will
            be a MultiIndex with `gene` and `taxon` as elements. If the equals
            'gene', the row index will be gene based and if 'taxon' will be
            taxon based
        gene_func (func): a function to map a gene_id to a gene_map. See
            :func:`.mapper.map_gene_id` for an example
        taxon_func (func): a function to map a taxon_id to a list of IDs. See
            :mod:`.mapper.map_taxon_id_to_rank` or
            :mod:`.mapper.map_taxon_id_to_ancestor` for examples

    Returns:
        DataFrame: :class:`pandas.DataFrame` with the pN/pS values for the input
        SNPs, with the columns being the samples.

    """
    sample_dict = dict((sample, {}) for sample in snps_data)
    multi_index = set()

    if gene_func is None:
        gene_func = functools.partial(itertools.repeat, times=1)
    if taxon_func is None:
        taxon_func = functools.partial(itertools.repeat, times=1)

    for sample, genes_dict in snps_data.iteritems():

        LOG.info('Analysing SNP from sample %s', sample)

        for gene_syn in pipe_filters(genes_dict.itervalues(), *filters):

            #in old data the gene_id was in the form K00001.UID
            #it's now allowed anymore
            if isinstance(gene_syn, GeneSyn):
                gene_syn.gene_id = gene_syn.gene_id.split('.')[0]

            iter_func = itertools.product(
                gene_func(gene_syn.gene_id),
                taxon_func(gene_syn.taxon_id)
            )

            for gene_id, taxon_id in iter_func:

                if index_type == 'gene':
                    key = gene_id
                elif index_type == 'taxon':
                    key = taxon_id
                else:
                    key = (gene_id, taxon_id)

                #we don't care about info about ids and so on, only syn/nonsyn
                #and coverage, to use the calc_ratio method
                try:
                    sample_dict[sample][key].add(gene_syn)
                except KeyError:
                    #Needed with the new GeneSNP, because copies the references
                    #and the number of SNPs raises (the original data structure
                    #is modified)
                    sample_dict[sample][key] = copy.deepcopy(gene_syn)

                multi_index.add(key)

    if isinstance(key, tuple):
        multi_index = pandas.MultiIndex.from_tuples(
            sorted(multi_index),
            names=('gene', 'taxon')
        )
    else:
        multi_index = pandas.Index(multi_index)

    #we already satisfied a minimum coverage filter or at least if doesn't
    #matter in the calculation anymore, using haplotypes=True, the special
    #case where syn=nonsyn=0 will result in a 0 as pN/pS for a GeneSyn instance
    sample_dict = dict(
        (
            sample,
            dict(
                (key, gene.calc_ratio(haplotypes=True))
                for key, gene in row_dict.iteritems()
            )
        )
        for sample, row_dict in sample_dict.iteritems()
    )
    dataframe = pandas.DataFrame(sample_dict, index=multi_index,
                                 columns=sorted(sample_dict.keys()))

    return dataframe[dataframe.count(axis=1) >= min_num]


def significance_test(dataframe, taxon_id1, taxon_id2,
                      test_func=scipy.stats.ks_2samp):
    """
    .. versionadded:: 0.1.11

    Perform a statistical test on each gene distribution in two different taxa.

    For each gene common to the two taxa, the distribution of values in all
    samples (columns) between the two specified taxa is tested.

    Arguments:
        dataframe: :class:`pandas.DataFrame` instance
        taxon_id1: the first taxon ID
        taxon_id2: the second taxon ID
        test_func: function used to test,
          defaults to :func:`scipy.stats.ks_2samp`

    Returns:
        pandas.Series: with all pvalues from the tests
    """
    LOG.info("Performing %s test", test_func.__name__)

    gene_ids = set(
        dataframe.select(
            lambda x: x[1] == taxon_id1
        ).index.get_level_values('gene') &
        dataframe.select(
            lambda x: x[1] == taxon_id2
        ).index.get_level_values('gene')
    )
    LOG.info("Number of genes to be tested: %d", len(gene_ids))

    sign_genes = {}

    for gene_id in gene_ids:
        tx1_vals = dataframe.loc[gene_id].loc[taxon_id1].dropna()
        tx2_vals = dataframe.loc[gene_id].loc[taxon_id2].dropna()

        pvalue = test_func(tx1_vals, tx2_vals)[1]

        sign_genes[gene_id] = pvalue

    return pandas.Series(sign_genes)


def flat_sample_snps(snps_data, min_cov):
    """
    .. versionadded:: 0.1.11

    Adds all the values of a gene across all samples into one instance of
    :class:`classes.GeneSyn`, giving the average gene among all samples.

    Arguments:
        snps_data (dict): dictionary with the instances of
            :class:`classes.GeneSyn`
        min_cov (int): minimum coverage required for the each instance to be
            added

    Returns:
        dict: the dictionary with only one key (`all_samples`), which can be
        used with :func:`combine_sample_snps`
    """
    samples = snps_data.keys()
    gene_ids = snps_data[samples[0]].keys()
    new_data = {'all_samples': {}}

    for gene_id in gene_ids:
        for sample in samples:
            gene_syn = snps_data[sample][gene_id]

            if gene_syn.coverage < min_cov:
                continue

            try:
                new_data['all_samples'][gene_id].add(gene_syn)
            except KeyError:
                new_data['all_samples'][gene_id] = copy.copy(gene_syn)

    return new_data
