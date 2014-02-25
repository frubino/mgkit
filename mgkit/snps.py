# encoding=utf8
"""
Manage SNP data.

"""
import logging
import itertools
import numpy
import pandas
import scipy.stats
import csv
from .taxon import BLACK_LIST

LOG = logging.getLogger(__name__)

MIN_COV = 4
"Minumum coverage required in some functions."

MIN_NUM = 10
"Used to set the minimum number of replicates for some functions"


class GeneSyn(object):
    """
    Class defining gene and synonymous/non-synonymous SNPs.

    It defines background synonymous/non-synonymous attributes and only has a
    method right now, which calculate pN/pS ratio.

    Attributes:
        gene_id (str): gene id
        taxon_id (int): gene taxon
        exp_syn (int): expected synonymous changes
        exp_nonsyn (int): expected non-synonymous changes
        syn (int): synonymous changes
        nonsyn (int): non-synonymous changes
        coverage (int): gene coverage

    .. warning::

        the `gid` and `taxon` attributes will be renamed in `gene_id` and
        `taxon_id` in later versions of the library, so they shouldn't be used.

    """
    __slots__ = (
        'gid',
        'taxon',
        'exp_syn',
        'exp_nonsyn',
        'syn',
        'nonsyn',
        'coverage',
        'taxon_root'
    )

    def __init__(self, gid='', taxon='', exp_syn=0, exp_nonsyn=0, syn=0,
                 nonsyn=0, coverage=None, taxon_root='', gene_id='',
                 taxon_id=0):

        self.gid = gid
        self.taxon = taxon
        self.taxon_root = taxon_root
        self.exp_syn = exp_syn
        self.exp_nonsyn = exp_nonsyn
        self.syn = syn
        self.nonsyn = nonsyn
        self.coverage = coverage
        self.gid = gene_id
        self.taxon = taxon_id

    def calc_ratio(self, flag_value=False, min_cov=None, haplotypes=False):
        """
        Calculate :math:`\\frac {pN}{pS}` for the gene.

        .. math::
            :label: pn-ps

            \\frac {pN}{pS} = \\frac{ ^{oN}/_{eN}}{ ^{oS}/_{eS}}

        WHere:

        * oN (number of non-synonymous - **nonsyn**)
        * eN (expected number of non-synonymous - **exp_nonsyn**)
        * oS (number of synonymous - **syn**)
        * eS (expected number of synonymous - **exp_syn**)

        Arguments:
            flag_value (bool): when there's no way to calculate the ratio, the
                possible cases will be flagged with a negative number. This
                allows to make substitutions for these values
            min_cov (int, None): minimum coverage require for some special
            cases. if is None, it's set to the global variable :data:`MIN_COV`.
            haplotypes (bool): if true, coverage information is not used,
                because the SNPs are assumed to come from an alignment that has
                sequences having haplotypes

        Returns:
            float: the :math:`\\frac {pN}{pS}` for the gene.

            .. note::

                Because pN or pS can be 0, and the return value would be NaN, we
                take in account some special cases. The default return value in
                this cases is :const:`numpy.nan`.

            * Both synonymous and non-synonymous values are 0:

                * if both the syn and nonsyn attributes are 0 but there's
                  coverage for this gene, we return a 0, as there's no evolution
                  in this gene. The attribute coverage must be equal or greater
                  than the min_cov parameter, which is by default assigned from
                  the :data:`MIN_COV` global variable (if left at *None*): this
                  means that it can be configured at runtime, if we want a
                  different default value without passing it to the method for
                  each call
                * In case the **coverage** attribute is **None** and the
                  **flag_value** parameter is True, the return value is **-3**

            * The number of non-synonymous is greater than 0 but the number of
              synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-1**

            * The number of synonymous is greater than 0 but the number of
              non-synonymous is 0:

                * if **flag_value** is **True**, the returned value is **-2**

            +------------+------------+------------+----------+--------------+
            | :math:`oS` | :math:`oN` | flag_value | coverage | return value |
            +============+============+============+==========+==============+
            | 0          | 0          | Not Used   | `int`    | **0**        |
            +------------+------------+------------+----------+--------------+
            | 0          | 0          | True       | Not Used | **-3**       |
            +------------+------------+------------+----------+--------------+
            | >0         | 0          | True       | Not Used | **-1**       |
            +------------+------------+------------+----------+--------------+
            | 0          | >0         | True       | Not Used | **-2**       |
            +------------+------------+------------+----------+--------------+

        """
        #set the minimum coverage value if not specified
        if min_cov is None:
            min_cov = MIN_COV

        #Both values are non-zero
        if (self.nonsyn != 0) and (self.syn != 0):
            pN = self.nonsyn / float(self.exp_nonsyn)
            pS = self.syn / float(self.exp_syn)
            return pN / pS
        #case in which a the SNPs come from haplotypes, in this case we don't
        #need to check for coverage to return a 0 for this special case
        elif (self.nonsyn == 0) and (self.syn == 0) and haplotypes:
            return 0

        #case in which a coverage attribute is specified, in this case we don't
        #need to flag the return value
        #getattr is used in case the value was not initialized, in cases like
        #a deserialized instance had no coverage value
        if getattr(self, 'coverage', None) is not None:
            if (self.nonsyn == 0) and (self.syn == 0):
                if self.coverage >= min_cov:
                    # LOG.debug("Coverage (%d) OK", min_cov)
                    return 0
                # LOG.debug("Coverage (%d) KO", min_cov)

        if flag_value:
            if self.nonsyn != 0:
                #there's at least non-synonymous count but no synonymous one
                #this will be converted in the max value for the matrix or
                #to some other value (-1 flag this case)
                if self.syn == 0:
                    return -1
            else:
                #there's at least a synonymous count but no non-synonymous one
                #this will be converted in the max value for the matrix or
                #to some other value (-2 flag this case)
                if self.syn != 0:
                    return -2
                else:
                    #There's no changes in the gene at at all. It should be
                    #checked if that gene has coverage.
                    return -3

        return numpy.nan

    def __getstate__(self):
        return dict((x, getattr(self, x)) for x in self.__slots__)

    def __setstate__(self, state):
        for name, value, in state.iteritems():
            setattr(self, name, value)

    def to_string(self):
        """
        Return a string with some info about the instance. Used by __str__
        """
        return '{0}-{1} pN/pS: {2:.2f}'.format(
            self.gid, self.taxon if self.taxon else None, self.calc_ratio()
        )

    @property
    def gene_id(self):
        "Alias for gid attribute at the moment"
        return self.gid

    @gene_id.setter
    def gene_id(self, gene_id):
        self.gid = gene_id

    @property
    def taxon_id(self):
        "Alias for taxon attribute at the moment"
        return self.taxon

    @taxon_id.setter
    def taxon_id(self, taxon_id):
        self.taxon = taxon_id

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string()

    def add(self, other):
        """
        Inplace addition of another instance values. No check for them being the
        same gene/taxon, it's up to the user to check that they can be added
        together.

        Arguments:
            other: instance of :class:`GeneSyn` to add
        """
        self.exp_nonsyn += other.exp_nonsyn
        self.exp_syn += other.exp_syn
        self.syn += other.syn
        self.nonsyn += other.nonsyn
        #only adds up coverage if the attribute is at least present in the other
        #object.
        if other.coverage is not None:
            if self.coverage is None:
                self.coverage = other.coverage
            else:
                self.coverage += other.coverage


def combine_snps_in_dataframe_test(count_dict, taxonomy, min_cov=MIN_COV,
                              black_list=None, min_num=MIN_NUM, rank=None,
                              anc_map=None, rooted=True, var_map=None,
                              feature='gene-taxon', only_rank=False,
                              gene_map=None, taxa_filter=None):
    """
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
        * refactor the two functions to make it smaller/easier to maintain
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
                    if not key in sample_dict:
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

                    if not key in sample_dict:
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


def combine_snps_in_dataframe(count_dict, taxonomy, min_cov=MIN_COV,
                              black_list=None, min_num=MIN_NUM, rank=None,
                              anc_map=None, rooted=True, var_map=None,
                              feature='gene-taxon', only_rank=False):
    """
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

            if rank is not None:
                gene_taxon = taxonomy.get_ranked_taxon(
                    gene_taxon,
                    rank
                ).taxon_id
                if only_rank:
                    if taxonomy[gene_taxon].rank != rank:
                        continue

            if anc_map is None:

                if feature == 'gene-taxon':
                    if rooted:
                        key = (gene_id, root_taxon, gene_taxon)
                    else:
                        key = (gene_id, gene_taxon)
                elif feature == 'taxon':
                    key = gene_taxon
                elif feature == 'gene':
                    key = gene_id

                multi_index.add(key)

                if not key in sample_dict:
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

                    if not key in sample_dict:
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
                if taxonomy.get_ranked_taxon(taxon_id, taxon_rank).rank == taxon_rank
            )
        )

    gene_index = sorted(set(dataframe.index.get_level_values('gene')))

    # print taxon_index
    rank_matrix = pandas.DataFrame(index=gene_index, columns=taxon_index)

    for gene_id in rank_matrix.index:
        gene_array = dataframe.loc[gene_id]
        # print gene_id, type(gene_array)
        for taxon_id, rank in zip(gene_array.index, scipy.stats.rankdata(gene_array)):
            if taxon_rank is not None:
                taxon_id = taxonomy.get_ranked_taxon(taxon_id, taxon_rank).taxon_id
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
    rank_matrix = pandas.DataFrame(index=gene_map.keys(), columns=dataframe.columns)

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


def wilcoxon_pairwise_test_dataframe(dataframe, test_func=scipy.stats.ranksums,
                                     threshold=0.01):
    """
    Make a wilcoxon pairwise test on a dataframe whose rows are a MultiIndex
    object, in the order gene, taxon.

    For each gene, each taxon variant is tested against each other, in a
    pairwise manner. If at least one of the tests is significant, the gene is
    reported as significantly different.

    :dataframe: dataframe instance
    :test_func: function used to test, defaults to :func:`scipy.stats.ranksums`
    :threshold: threshold for null hypotesis discard.

    .. note::

        The function access the row using the
        'loc' attribute on the object, so the order must be correct. No checks
        are made at the moment on the object for the consistence of it.
    """
    LOG.info("Performing Wilcoxon ranksums test pairwise")

    genes = set(dataframe.index.get_level_values('gene'))
    LOG.info("Number of genes to be tested: %d", len(genes))

    sign_genes = []

    for gene in genes:
        dataframe_gene = dataframe.loc[gene]

        if len(dataframe_gene) <= 1:
            continue

        for taxon_id1, taxon_id2 in itertools.combinations(dataframe_gene.index, 2):
            pvalue = test_func(dataframe_gene.loc[taxon_id1].dropna(),
                               dataframe_gene.loc[taxon_id2].dropna())[1]

            if pvalue < threshold:
                sign_genes.append(gene)
                break

    return sign_genes


def write_sign_genes_table(out_file, dataframe, sign_genes, taxonomy,
                           gene_names=None, sep=','):
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


import copy

from .filter.snps import pipe_filters


def combine_sample_snps(snps_data, min_num, key_attrs, filters, gene_map=None,
                        taxon_map=None):
    """
    filtri da concaternare:
        * filters:
            * default:
                * coverage
                * black list taxa
            * taxa_filter
        * mapper:
            * gene_map
            * rank
            * anc_map


    """
    sample_dict = dict((sample, {}) for sample in snps_data)
    multi_index = set()

    for sample, genes_dict in snps_data.iteritems():
        print sample
        for gene_syn in pipe_filters(genes_dict.itervalues(), *filters):

            gene_syn.gene_id = gene_syn.gene_id.split('.')[0]

            key = tuple(
                getattr(gene_syn, attr)
                for attr in key_attrs
            )
            #no use in keeping a tuple if only one-element key
            if len(key) == 1:
                key = key[0]

            #we don't care about info about ids and so on, only syn/nonsyn
            #and coverage, to use the calc_ratio method
            try:
                sample_dict[sample][key].add(gene_syn)
            except KeyError:
                sample_dict[sample][key] = copy.copy(gene_syn)

            multi_index.add(key)

    if len(key_attrs) > 1:
        multi_index = pandas.MultiIndex.from_tuples(
            sorted(multi_index),
            names=tuple(attr[:-3] for attr in key_attrs)
        )
    else:
        multi_index = pandas.Index(multi_index)

    sample_dict = dict(
        (
            sample,
            dict(
                (key, gene.calc_ratio())
                for key, gene in row_dict.iteritems()
            )
        )
        for sample, row_dict in sample_dict.iteritems()
    )
    dataframe = pandas.DataFrame(sample_dict, index=multi_index,
                                 columns=sorted(sample_dict.keys()))

    return dataframe[dataframe.count(axis=1) >= min_num]
