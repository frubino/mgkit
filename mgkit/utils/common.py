"""
Utility functions
"""
from __future__ import division
import itertools
import operator
import functools
import warnings


def average_length(a1s, a1e, a2s, a2e):
    """
    Given two sets of coordinates, a1 and a2, returns the average length.

    :param int a1s: a1 leftmost number
    :param int a1e: a1 rightmost number
    :param int a2s: a2 leftmost number
    :param int a2e: a2 rightmost number

    :return float: the average length
    """
    return ((a1e - a1s) + (a2e - a2s) + 2) / 2.0


def between(pos, start, end):
    """
    Tests if a number is between two others

    :param int pos: number to test
    :param int start: leftmost number
    :param int end: rightmost number

    :return bool: if the number is between start and end
    """
    if pos < start or pos > end:
        return False
    return True


def union_range(start1, end1, start2, end2):
    """
    .. versionadded:: 0.1.12

    If two numeric ranges overlap, it returns the new range, otherwise None is
    returned. Works on both int and float numbers, even mixed.

    Arguments:
        start1 (numeric): start of range 1
        end1 (numeric): end of range 1
        start2 (numeric): start of range 2
        end2 (numeric): end of range 2

    Returns:
        (tuple or None): union of the ranges or None if the ranges don't overlap

    Example:
        >>> union_range(10, 13, 1, 10)
        (1, 13)
        >>> union_range(1, 10, 11, 13)
        None

    """
    if between(start2, start1, end1) or between(end2, start1, end1):
        return min(start1, start2), max(end1, end2)
    return None


def ranges_length(ranges):
    """
    .. versionadded:: 0.1.12

    Given an iterable where each element is a range, a tuple whose elements
    are numbers with the first being less than or equal to the second, the
    function sums the lengths of all ranges.

    Arguments:
        ranges (iterable): each element is a tuple like *(1, 10)*

    Returns:
        int: sum of all ranges lengths

    """
    return sum(range[1] - range[0] + 1 for range in ranges)


def range_substract(start1, end1, start2, end2):
    intersect = range_intersect(start1, end1, start2, end2)

    if intersect is None:
        return [(start1, end1)]

    ranges = []
    if start1 != intersect[0]:
        ranges.append((start1, intersect[0] - 1))
    if end1 != intersect[1]:
        ranges.append((intersect[1] + 1, end1))
    return ranges


def range_intersect(start1, end1, start2, end2):
    """
    .. versionadded:: 0.1.13

    Given two ranges in the form *(start, end)*, it returns the range
    that is the intersection of the two.

    Arguments:
        start1 (int): start position for the first range
        end1 (int): end position for the first range
        start2 (int): start position for the second range
        end2 (int): end position for the second range

    Returns:
        (None, tuple): returns a tuple with the start and end position for
        the intersection of the two ranges, or *None* if the intersection is
        empty
    """
    if between(start2, start1, end1) or between(end2, start1, end1) or \
       between(start1, start2, end2) or between(end1, start2, end2):
        return max(start1, start2), min(end1, end2)
    return None


def apply_func_window(func, data, window, step=0):

    if step == 0:
        step = window

    for index in xrange(0, len(data), step):
        yield func(data[index:index+window])


def range_substract_(range1, ranges):
    range1 = set(xrange(range1[0], range1[1] + 1))
    range1.difference_update(
        *(
            xrange(x[0], x[1] + 1)
            for x in ranges
        )
    )

    for k, g in itertools.groupby(enumerate(range1), lambda (i, x): i - x):
        group = map(operator.itemgetter(1), g)
        yield min(group), max(group)


def deprecated(func):
    '''
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.

    from https://wiki.python.org/moin/PythonDecoratorLibrary
    '''

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
            "Call to deprecated function {}.\n{}".format(
                func.__name__,
                func.__doc__
            ),
            category=DeprecationWarning,
            filename=func.func_code.co_filename,
            lineno=func.func_code.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func
