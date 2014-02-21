"""
Module used to map taxon_id to different levels in the taxonomy.
"""


def map_taxon_by_id_list(taxon_id, map_ids, func):
    """
    Maps a taxon_id to a list of taxon IDs, using the function supplied.

    Arguments:
        taxon_id (int): taxon ID to map
        map_ids (iterable): list of taxon IDs to which the taxon_id will be
            mapped.
        func (func): function used to map the IDs, accepts two taxon IDs

    Results:
        generator: generator expression of all IDs in map_ids to which taxon_id
            can be mapped.

    Example:
        If mapping a taxon (Prevotella ruminicola) to Prevotella or Clostridium,
        using as `func` :func:`mgkit.taxon.is_ancestor` and taxonomy is an
        instance of :class:`mgkit.taxon.UniprotTaxonomy`.

        >>> import functools
        >>> from mgkit.taxon import is_ancestor
        >>> func = functools.partial(is_ancestor, taxonomy)
        >>> list(map_taxon_by_id_list(839, [838, 1485], func))
        [838]

    """
    taxon_ids = (
        map_id
        for map_id in map_ids
        if func(taxon_id, map_id)
    )

    return taxon_ids
