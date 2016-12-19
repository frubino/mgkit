"""
Various utilities to help read and process files
"""
import logging


LOG = logging.getLogger(__name__)


def group_tuples_by_key(iterator, key_func=None, skip_elements=0):
    """
    .. versionadded:: 0.3.1

    Group the elements of an iterator by a key and yields the grouped elements.
    The elements yielded by the iterator are assumed to be a list or tuple,
    with the default key (when *key_func* is None) being the first of the of
    the objects inside that element. This behaviour can be customised by
    passing to *key_func* a function that accept an element and returns the key
    to be used.

    .. note::

        the iterable assumen that the elements are already sorted by their keys

    Arguments:
        iterator (iterable): iterator to be grouped
        key_func (func): function that accepts a element and returns its
            associated key
        skip_elements (int): number of elements to skip at the start

    Yields:
        list: a list of the grouped elements by key
    """
    if key_func is None:
        key_func = lambda x: x[0]

    for index in xrange(skip_elements):
        next(iterator)

    curr_key = None
    curr_ann = []

    for element in iterator:
        new_key = key_func(element)
        if curr_key == new_key:
            curr_ann.append(element)
        else:
            if curr_key is None:
                curr_ann.append(element)
                curr_key = new_key
            else:
                yield curr_ann
                curr_key = new_key
                curr_ann = [element]
    else:
        yield curr_ann
