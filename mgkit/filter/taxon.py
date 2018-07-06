"""
.. versionadded:: 0.1.9

Taxa filtering functions
"""

from .common import FilterFails
from ..taxon import is_ancestor
import functools


def filter_taxon_by_id_list(taxon_id, filter_list=None, exclude=False,
                            func=None):
    """
    Filter a taxon_id against a list of taxon ids. Returns True if the
    conditions of the filter are met.

    If func is not None, a function that accepts two values is expected,
    it should be either a partial `is_ancestor` which only accepts taxon_id and
    anc_id or another function that behaves the same way.

    .. note::

            if func is None, a simple lambda is used to test identity::

                func = lambda t_id, a_id: t_id == a_id

    Arguments:
        taxon_id (int): the taxon id to filter
        filter_list (iterable): an iterable with taxon ids
        exclude (bool): if the filter is reversed (i.e. included if NOT found)
        func (func or None): a function that accepts taxon_id and an anc_id
            and returns a bool to indicated if anc_id is ancestor of taxon_id.
            Equivalent to :func:`~mgkit.taxon.is_ancestor`.

    Returns:
        bool:
            True if the taxon_id is in the filter list (or a descendant of it)
            False if it's not found. Exclude equal to True reverse the result.

            +-------+---------+--------------+
            | Found | Exclude | Return Value |
            +=======+=========+==============+
            | Yes   | False   | True         |
            +-------+---------+--------------+
            | No    | False   | False        |
            +-------+---------+--------------+
            | Yes   | True    | False        |
            +-------+---------+--------------+
            | No    | True    | True         |
            +-------+---------+--------------+

    Example:

        If using func and assuming that `taxonomy` is an instance of
        :class:`~mgkit.taxon.Taxonomy` with data loaded:

        >>> import functools
        >>> import mgkit.taxon
        >>> func = functools.partial(mgkit.taxon.is_ancestor, taxonomy)
        >>> filter_taxon_by_id_list(1200582, [838], func=func)
        True

    """

    if func is None:
        def func(t_id, a_id): return t_id == a_id

    if filter_list is None:
        raise FilterFails('No filter_list')

    ret_val = any(
        func(taxon_id, filter_id)
        for filter_id in filter_list
    )

    return ret_val ^ exclude


def filter_by_ancestor(taxon_id, filter_list=None, exclude=False,
                       taxonomy=None):
    """
    .. versionadded:: 0.1.13

    Convenience function for :func:`filter_taxon_by_id_list`, as explained in
    the latter example.
    """
    if filter_list is None:
        raise FilterFails('No filter_list')

    func = functools.partial(is_ancestor, taxonomy)

    return filter_taxon_by_id_list(
        taxon_id,
        filter_list=filter_list,
        func=func,
        exclude=exclude
    )
