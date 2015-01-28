"""
Wappers to use some of the general function of the snps package
in a simpler way.
"""
import functools
import mgkit.snps.funcs
import mgkit.snps.filter
import mgkit.snps.mapper


def get_rank_dataframe(snp_data, taxonomy, min_num=3, rank='order',
                       index_type='taxon'):
    """
    .. versionadded:: 0.1.11

    Returns a :class:`~pandas.DataFrame` with the pN/pS of the given
    SNPs data, mapping all taxa to the specified rank. Higher taxa won't
    be included.

    Shortcut for using :func:`~mgkit.snps.funcs.combine_sample_snps`, using
    filters from :func:`~mgkit.snps.filter.get_default_filters` and as
    `taxon_func` parameter :func:`~mgkit.snps.mapper.map_taxon_id_to_rank`,
    with include_higher equals to False

    Arguments:
        snp_data (dict): dictionary sample->GeneSyn of SNPs data
        taxonomy: Uniprot Taxonomy
        min_num (int): minimum number of samples in which a valid pN/pS
            is found
        rank (str): taxon rank to map. Valid ranks are found in
            :data:`mgkit.taxon.TAXON_RANKS`
        index_type (str, None): type of index to return

    Returns:
        DataFrame: :class:`pandas.DataFrame` of pN/pS values. The index type
        is 'taxon'
    """
    taxon_func = functools.partial(
        mgkit.snps.mapper.map_taxon_id_to_rank,
        taxonomy=taxonomy,
        rank=rank
    )

    filters = mgkit.snps.filter.get_default_filters(taxonomy)

    dataframe = mgkit.snps.funcs.combine_sample_snps(
        snp_data,
        min_num,
        filters,
        taxon_func=taxon_func,
        gene_func=None,
        index_type=index_type
    )

    return dataframe


def get_gene_map_dataframe(snp_data, taxonomy, gene_map, min_num=3,
                           index_type='gene'):
    """
    .. versionadded:: 0.1.11

    Returns a :class:`~pandas.DataFrame` with the pN/pS of the given
    SNPs data, mapping all taxa to the gene map.

    Shortcut for using :func:`~mgkit.snps.funcs.combine_sample_snps`, using
    filters from :func:`~mgkit.snps.filter.get_default_filters` and as
    `gene_func` parameter :func:`~mgkit.snps.mapper.map_gene_id`.

    Arguments:
        snp_data (dict): dictionary sample->GeneSyn of SNPs data
        taxonomy: Uniprot Taxonomy
        min_num (int): minimum number of samples in which a valid pN/pS
            is found
        gene_map (dict): dictionary of mapping for the gene_ids in in SNPs
            data
        index_type (str, None): type of index to return

    Returns:
        DataFrame: :class:`pandas.DataFrame` of pN/pS values. The index type
        is 'gene'
    """
    gene_func = functools.partial(
        mgkit.snps.mapper.map_gene_id,
        gene_map=gene_map
    )

    filters = mgkit.snps.filter.get_default_filters(taxonomy)

    dataframe = mgkit.snps.funcs.combine_sample_snps(
        snp_data,
        min_num,
        filters,
        taxon_func=None,
        gene_func=gene_func,
        index_type=index_type
    )

    return dataframe


def get_full_dataframe(snp_data, taxonomy, min_num=3, index_type=None):
    """
    .. versionadded:: 0.1.12

    Returns a :class:`~pandas.DataFrame` with the pN/pS of the given
    SNPs data.

    Shortcut for using :func:`~mgkit.snps.funcs.combine_sample_snps`, using
    filters from :func:`~mgkit.snps.filter.get_default_filters`.

    Arguments:
        snp_data (dict): dictionary sample->GeneSyn of SNPs data
        taxonomy: Uniprot Taxonomy
        min_num (int): minimum number of samples in which a valid pN/pS
            is found
        index_type (str, None): type of index to return

    Returns:
        DataFrame: :class:`pandas.DataFrame` of pN/pS values. The index type
        is None (gene-taxon)
    """

    filters = mgkit.snps.filter.get_default_filters(taxonomy)

    dataframe = mgkit.snps.funcs.combine_sample_snps(
        snp_data,
        min_num,
        filters,
        taxon_func=None,
        gene_func=None,
        index_type=index_type
    )

    return dataframe


def get_gene_taxon_dataframe(snp_data, taxonomy, gene_map, min_num=3,
                             rank='genus', index_type=None):
    """
    .. versionadded:: 0.1.12

    .. todo::

        edit docstring

    Returns a :class:`~pandas.DataFrame` with the pN/pS of the given
    SNPs data, mapping all taxa to the gene map.

    Shortcut for using :func:`~mgkit.snps.funcs.combine_sample_snps`, using
    filters from :func:`~mgkit.snps.filter.get_default_filters` and as
    `gene_func` parameter :func:`~mgkit.snps.mapper.map_gene_id`.

    Arguments:
        snp_data (dict): dictionary sample->GeneSyn of SNPs data
        taxonomy: Uniprot Taxonomy
        min_num (int): minimum number of samples in which a valid pN/pS
            is found
        gene_map (dict): dictionary of mapping for the gene_ids in in SNPs
            data
        index_type (str, None): type of index to return

    Returns:
        DataFrame: :class:`pandas.DataFrame` of pN/pS values. The index type
        is 'gene'
    """
    gene_func = functools.partial(
        mgkit.snps.mapper.map_gene_id,
        gene_map=gene_map
    )

    if rank is None:
        taxon_func = None
    else:
        taxon_func = functools.partial(
            mgkit.snps.mapper.map_taxon_id_to_rank,
            taxonomy=taxonomy,
            rank=rank
        )

    filters = mgkit.snps.filter.get_default_filters(taxonomy)

    dataframe = mgkit.snps.funcs.combine_sample_snps(
        snp_data,
        min_num,
        filters,
        taxon_func=taxon_func,
        gene_func=gene_func,
        index_type=index_type
    )

    return dataframe
