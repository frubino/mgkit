"""
Mapping functions for SNPs - Should be move into an 'iterator' package to
be shared with other modules?
"""


def map_gene_id(gene_id, gene_map=None):
    """
    Returns an iterator for all the values of a dictionary. if gene_id is not
    found in the gene_map, an empty iterator is returned.

    Arguments:
        gene_id (immutable): gene_id or any other dictionary key.
        gene_map (dict): a dictionary in the form key->[v1, v2, .. vN]

    Returns:
        generator: iterator (empty if gene_id is not in gene_map) with the
        values

    """
    if gene_id not in gene_map:
        return
    for map_id in gene_map[gene_id]:
        yield map_id


def map_taxon_id_to_rank(taxon_id, rank=None, taxonomy=None,
                         include_higher=False):
    """
    Given a taxon_id, returns an iterator with only the element that correspond
    to the requested rank. If the taxon returned by
    :class:`mgkit.taxon.UniprotTaxonomy.get_ranked_taxon` has a different rank
    than requested, the iterator will be empty if `include_higher` is False
    and the returned taxon ID if True.

    Arguments:
        taxon_id (int): taxon ID to be mapped
        rank (str): taxon rank used (:data:`mgkit.taxon.TAXON_RANKS`)
        include_higher (bool): determines if a rank higher than the one
            requested is to be returned

    Returns:
        generator: iterator with the values or empty
    """
    ranked_taxon = taxonomy.get_ranked_taxon(taxon_id, rank=rank)

    #in case we don't want to include higher taxa
    if (ranked_taxon.rank != rank) and (not include_higher):
        return

    yield ranked_taxon.taxon_id


def map_taxon_id_to_ancestor(taxon_id, anc_ids=None, func=None):
    """
    Given a taxon_id and a list of ancestors IDs, returns an iterator with the
    IDs that are ancestors of taxon_id.

    Arguments:
        taxon_id (int): taxon ID to be mapped
        anc_ids (iterable): taxon IDs to check for ancestry
        func: function used to check for ancestry - partial function for
            :func:`mgkit.taxon.is_ancestor` that accepts taxon_id and anc_id

    Returns:
        generator: iterator with the values or empty

    .. note::

        check :func:`mgkit.filter.taxon.filter_taxon_by_id_list` for examples
        on using func

    """
    for anc_id in anc_ids:
        if func(taxon_id, anc_id):
            yield anc_id
