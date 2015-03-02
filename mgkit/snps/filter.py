"""
SNPs filtering functions
"""
import functools
import itertools
from .. import consts
from ..filter.taxon import filter_taxon_by_id_list
from ..filter.common import FilterFails


def filter_genesyn_by_taxon_id(gene_syn, taxonomy=None, filter_list=None,
                               exclude=False, func=None):
    """
    Checks if the `taxon_id` attribute of `gene_syn` is the `filter_list`.
    Excelude reverses the result. If func is supplied, it's used to traverse the
    `taxonomy`.

    Arguments:
        gene_syn: :class:`~mgkit.snps.GeneSyn` instance
        taxonomy: a valid taxonomy (instance of
            :class:`~mgkit.taxon.UniprotTaxonomy`)
        filter_list (iterable): list of taxon IDs to include/exclude
        exclude (bool): if the filter is reversed
        func (func): :func:`~mgkit.taxon.is_ancestor`

    Returns:
        bool: if the exclude is True, the gene_id must appear in the gene_ids,
        if False, returns True only if gene_id is NOT in gene_ids.

    Raises:
        FilterFails: if filter_list is None or taxonomy is None and func is not
            None

    """

    if ((taxonomy is None) and (func is not None)) or (filter_list is None):
        raise FilterFails('No taxonomy supllied')

    if func is not None:
        func = functools.partial(func, taxonomy)

    taxon_id = gene_syn.taxon_id

    return filter_taxon_by_id_list(
        taxon_id,
        filter_list=filter_list,
        exclude=exclude,
        func=func
    )


def filter_genesyn_by_gene_id(gene_syn, gene_ids=None, exclude=False,
                              id_func=None):
    """
    Checks if the gene_id is listed in the filter_list.

    Arguments:
        gene_syn: :class:`~mgkit.snps.GeneSyn` instance
        gene_ids (iterable): list of gene IDs to include/exclude
        exclude (bool): if the filter is reversed

    Returns:
        bool: if the exclude is True, the gene_id must appear in the gene_ids,
        if False, returns True only if gene_id is NOT in gene_ids.

    Raises:
        FilterFails: if gene_ids is None

    """
    if gene_ids is None:
        raise FilterFails('No gene_ids supplied')
    return (id_func(gene_syn) in gene_ids) ^ exclude


def filter_genesyn_by_coverage(gene_syn, min_cov=None):
    """
    Checks if the coverage of the provided `gene_syn` is at least `min_cov`

    Arguments:
        gene_syn: :class:`~mgkit.snps.GeneSyn` instance
        min_cov (int): minimum coverage allowed (included)

    Returns:
        bool: True if the gene has enough coverage

    Raises:
        FilterFails: if min_cov is None
    """
    if min_cov is None:
        raise FilterFails('No coverage supplied')

    return gene_syn.coverage >= min_cov


def get_default_filters(taxonomy, **kwargs):
    """
    Retuns a list of filters that are used by default. it needs a valid taxonomy
    and gets the default arguments from :data:`mgkit.consts.DEFAULT_SNP_FILTER`.
    """
    filter_opts = consts.DEFAULT_SNP_FILTER.copy()
    filter_opts.update(kwargs)

    filter_coverage = functools.partial(
        filter_genesyn_by_coverage,
        min_cov=filter_opts['min_cov']
    )
    filter_black_list = functools.partial(
        filter_genesyn_by_taxon_id,
        taxonomy=taxonomy,
        filter_list=filter_opts['include_only'],
        exclude=False,
        func=filter_opts['func']
    )
    return [filter_coverage, filter_black_list]


def pipe_filters(iterable, *funcs):
    """
    Pipes a list of filter to iterable, using the python ifilter function in
    the itertools module.
    """
    for func in funcs:
        iterable = itertools.ifilter(func, iterable)
    for value in iterable:
        yield value
