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
from ..consts import BLACK_LIST
from .classes import GeneSyn
from .. import consts
from .filter import pipe_filters

LOG = logging.getLogger(__name__)


def combine_snps_in_dataframe_test(count_dict, taxonomy, min_cov=consts.MIN_COV,
                                   black_list=None, min_num=consts.MIN_NUM,
                                   rank=None, anc_map=None, rooted=True,
                                   var_map=None, feature='gene-taxon',
                                   only_rank=False, gene_map=None, taxa_filter=None):
    """
    .. deprecated:: 0.1.11
        use :func:`combine_sample_snps` instead

    Convert a sample->gene->GeneSyn dictionary into a :class:`pandas.DataFrame`
    with a :class:`pandas.MultiIndex` composed of gene_id, root_taxon, taxon as
    row index and an :class:`pandas.Index` with the sample names as column
    index.

    :param dict count_dict: dictionary containing :class:`GeneSyn` instances of
        KO_IDX sample->ko_idx->:class:`GeneSyn`
    :param dict tmap: taxon root map obtained from :func:`taxon.group_by_root`
        with only_names=True. If None the root won't be added to the index
    :param int min_cov: minimum coverage allowed for inclusion of a gene
    :param list black_list: black list of taxa to exclude; defaults to
        :data:`taxon.BLACK_LIST`
    :param int min_num: minimum number of replicates allowed for inclusion in
        the resulting matrix
    :param str rank: taxon level that represent the maximum specifity of the
        dataset
    :param bool rooted: if True, the index will include the root taxon of the
        gene-taxon tuple (as a string)
    :param dict anc_map: dictionary with the ancestry map to match a profile
        based taxonomy
    :param dict var_map: dictionary instance, used only if anc_map is not None,
        to store the ko_idx that belong to a profile
    :param str feature: one of 'gene', 'taxon', 'gene-taxon', works in the same
        way as combine_genes_by_feature
    :param bool only_rank: only active if rank is specified. If True only the
        defined rank is included in the results. Any rank above is taken out.

    .. note::

        if anc_map is defined, the number of rows in the resulting DataFrame
        may not correspond to the actual number of profiles, because the
        DataFrame building doesn't take in account the genes that are in the
        profiles (the guide is just the taxa), anc_map is returned by
        :func:`taxon.get_ancestor_map`

    :return DataFrame: returns a :class:`pandas.DataFrame`

    .. todo::

        * eliminate this duplicate function
    """
    LOG.info("Combining gene variants into DataFrame")
    LOG.info("Feature selected: %s", feature)
    if black_list is None:
        black_list = BLACK_LIST
    LOG.info("Black listed taxa: %s", ','.join(black_list))

    if anc_map is not None:
        LOG.info("Using ancestry map")
    if not rooted:
        LOG.info("Using index with no root")

    sample_dict = dict((sample, {}) for sample in count_dict)
    multi_index = set()

    skipped = 0
    added = 0

    for sample, genes_dict in count_dict.iteritems():
        LOG.debug("Sample %s loop", sample)
        for gene in genes_dict.itervalues():

            gene_taxon = gene.taxon

            #skips genes belonging to taxa black listed
            if any(True for blisted in black_list
                   if blisted in taxonomy[gene_taxon].lineage or
                   taxonomy[gene_taxon].s_name == blisted):
                # print "Skipped", taxonomy[gene_taxon].s_name
                continue

            if rooted:
                root_taxon = taxonomy.get_taxon_root(gene.taxon).s_name
            else:
                root_taxon = None

            #skips genes that have low coverage
            if gene.coverage < min_cov:
                skipped += 1
                continue
            added += 1

            gene_id = gene.gid.split('.')[0]

            #make the id an iterable, to conform to what mappings usually are
            if gene_map is None:
                gene_ids = [gene_id]
            else:
                #if the gene can't be mapped is skipped
                try:
                    gene_ids = gene_map[gene_id]
                except KeyError:
                    continue

            if rank is not None:
                gene_taxon = taxonomy.get_ranked_taxon(
                    gene_taxon,
                    rank
                ).taxon_id
                if only_rank:
                    if taxonomy[gene_taxon].rank != rank:
                        continue

            if taxa_filter is not None:
                if gene_taxon not in taxa_filter:
                    continue
                # LOG.debug("Inluding taxa %r", taxonomy[gene_taxon])

            if anc_map is None:

                if feature == 'gene-taxon':
                    if rooted:
                        keys = [
                            (gene_id, root_taxon, gene_taxon)
                            for gene_id in gene_ids
                        ]
                    else:
                        keys = [
                            (gene_id, gene_taxon)
                            for gene_id in gene_ids
                        ]
                elif feature == 'taxon':
                    keys = [gene_taxon]
                elif feature == 'gene':
                    keys = gene_ids

                multi_index.update(keys)

                for key in keys:
                    if not key in sample_dict[sample]:
                        sample_dict[sample][key] = GeneSyn(
                            gid=key,
                            taxon=gene_taxon,
                            taxon_root=root_taxon
                        )

                    sample_dict[sample][key].add(gene)

            else:
                #when the rows that have at least min_num of values are dropped
                #the index is unchanged. So some can be still in the anc_map,
                #which is based off the index.
                if gene_taxon not in anc_map:
                    continue

                for anc_taxon in anc_map[gene_taxon]:

                    if var_map is not None:
                        var_key = (gene_id, anc_taxon)
                        try:
                            var_map[var_key].add(gene.gid)
                        except KeyError:
                            var_map[var_key] = set()
                            var_map[var_key].add(gene.gid)

                    if rooted:
                        key = (gene_id, root_taxon, anc_taxon)
                    else:
                        key = (gene_id, anc_taxon)

                    multi_index.add(key)

                    if not key in sample_dict[sample]:
                        sample_dict[sample][key] = GeneSyn(
                            gid=key,
                            taxon=anc_taxon,
                            taxon_root=root_taxon
                        )

                    sample_dict[sample][key].add(gene)

    if feature == 'gene-taxon':
        multi_index = pandas.MultiIndex.from_tuples(
            sorted(multi_index),
            names=('gene', 'root', 'taxon') if rooted else ('gene', 'taxon')
        )
    else:
        multi_index = pandas.Index(multi_index)

    LOG.debug("Calculating pN/pS for all matrix elements")

    sample_dict = dict(
        (
            sample,
            dict(
                (idx, gene.calc_ratio())
                for idx, gene in row_dict.iteritems()
            )
        )
        for sample, row_dict in sample_dict.iteritems()
    )

    LOG.debug("Skipped genes for coverage (%d): %.2f%%", min_cov,
              float(skipped) / (skipped + added) * 100)

    dataframe = pandas.DataFrame(sample_dict, index=multi_index,
                                 columns=sorted(sample_dict.keys()))

    return dataframe[dataframe.count(axis=1) >= min_num]


def get_values_partition(observed, profile, neutral=1.0):
    "Return the partition to which an observed value belongs"

    in_c_or_d = lambda obs, prof: (obs < neutral) and (prof < neutral)

    funcs = {
        'A': lambda obs, prof: (obs > neutral) and (prof < neutral),
        'B': lambda obs, prof: (obs > neutral) and (prof > neutral),
        'C': lambda obs, prof: in_c_or_d(obs, prof) and (obs / prof > 1.0),
        'D': lambda obs, prof: in_c_or_d(obs, prof) and (obs / prof < 1.0),
        'E': lambda obs, prof: (obs < neutral) and (prof > neutral)
    }

    for partition, func in funcs.iteritems():
        if func(observed, profile):
            return partition

    return 'N'


def group_pnps_values_dataframe(observed_data, profile_data, taxonomy):
    """
    :param Series rumen_dict: a :class:`pandas.Series` with the index
        in the form (gene, taxon, root) and as values the mean pN/pS across all
        samples
    :param dict profile_dict: a dictionary with all the pN/pS values of the
        profiles with the key in the form (gene, taxon)
    :param iterable roots: list of root taxa to include, if None, all roots in
        the DataFramce are used

    :return dict: dictionary root->DataFrame
    """

    roots = set(
        taxonomy.get_taxon_root(taxon_id).s_name
        for taxon_id in set(observed_data.index.get_level_values('taxon'))
    )

    root_dict = dict((root, {'observed': {}, 'profile': {}}) for root in roots)

    for (gene_id, taxon_id), profile_value in profile_data.iterkv():
        gt_key = (gene_id, taxon_id)
        if gt_key not in observed_data:
            continue

        observed_value = observed_data[gt_key]
        #skips genes for which we have a NaN in the profile
        if numpy.isnan(profile_value) or numpy.isnan(observed_value):
            continue

        root = taxonomy.get_taxon_root(taxon_id).s_name

        root_dict[root]['observed'][gt_key] = observed_value
        root_dict[root]['profile'][gt_key] = profile_value

    root_dict = dict(
        (
            root, pandas.DataFrame(
                values,
                index=sorted(values['observed'].keys()),
                columns=('profile', 'observed')
            )
        )
        for root, values in root_dict.iteritems()
    )

    return root_dict


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
        for taxon_id, rank in mapped_matrix.mean().dropna().iterkv():
            rank_matrix.set_value(mapping_id, taxon_id, rank)

    return rank_matrix


def write_pnps_grouped_by_pathway(data, gene_list, taxa_list=None):
    """
    Groups a DataFrame values, where the index is a tuple (ko_id, taxon_id) and
    the the columns are profile and observed. The function group only the gene
    list provided (and optionally the taxa) into the partitions returned by
    :func:`get_values_partition`.
    """
    part_data = {}
    for (ko_id, taxon_id), (p_val, o_val) in data.iterrows():
        if (taxa_list is not None) and (taxon_id not in taxa_list):
            continue
        if ko_id in gene_list:
            part = get_values_partition(o_val, p_val)
            try:
                part_data[part].add(ko_id)
            except KeyError:
                part_data[part] = set([ko_id])

    return part_data


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
        snps_data (dict): dictionary with the `GeneSyn` instances
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
        print sample
        for gene_syn in pipe_filters(genes_dict.itervalues(), *filters):

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
                    sample_dict[sample][key] = copy.copy(gene_syn)

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
