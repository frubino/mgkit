"""
Utility functions
"""
from __future__ import division
from builtins import range
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

    .. versionchanged:: 0.3.1
        changed behaviour, since the intervals are meant to be closed

    If two numeric ranges overlap, it returns the new range, otherwise None is
    returned. Works on both int and float numbers, even mixed.

    Arguments:
        start1 (numeric): start of range 1
        end1 (numeric): end of range 1
        start2 (numeric): start of range 2
        end2 (numeric): end of range 2

    Returns:
        (tuple or None): union of the ranges or None if the ranges don't
        overlap

    Example:
        >>> union_range(10, 13, 1, 10)
        (1, 13)
        >>> union_range(1, 10, 11, 13)
        (1, 13)
        >>> union_range(1, 10, 12, 13)
        None

    """
    if between(start2, start1, end1) or \
        between(end2, start1, end1) or \
            (end1 + 1 == start2) or (end2 + 1 == start1):
        return min(start1, start2), max(end1, end2)
    return None


def union_ranges(intervals):
    """
    .. versionadded:: 0.3.1

    From a list of ranges, assumed to be closed, performs a union of all
    elements.

    Arguments:
        intervals (intervals): iterable where each element is a closed range
            (tuple)

    Returns:
        list: the list of ranges that are the union of all elements passed

    Examples:
        >>> union_ranges([(1, 2), (3, 7), (6, 12), (9, 17), (18, 20)])
        [(1, 20)]
        >>> union_ranges([(1, 2), (3, 7), (6, 12), (9, 14), (18, 20)])
        [(1, 14), (18, 20)]
    """
    intervals = sorted(intervals)

    union = [intervals.pop(0)]

    for start2, end2 in intervals:
        start1, end1 = union[-1]
        new_range = union_range(start1, end1, start2, end2)
        if new_range is None:
            union.append(
                (start2, end2)
            )
        else:
            union[-1] = new_range
    return union


def complement_ranges(intervals, end=None):
    """
    .. versionadded:: 0.3.1

    Perform a complement operation of the list of intervals, i.e. returning the
    ranges (tuples) that are not included in the list of intervals.
    :func:`union_ranges` is first called on the intervals.

    .. note::

        the `end` parameter is there for cases where the ranges passed don't
        cover the whole space. Assuming a list of ranges from annotations on a
        nucleotidic sequence, if the last range doesn't include the last
        position of the sequence, passing end equal to the length of the
        sequence will make the function include a last range that includes it

    Arguments:
        intervals (intervals): iterable where each element is a closed range
            (tuple)
        end (int): if the end of the complement intervals is supposed to be
            outside the last range.

    Returns:
        list: the list of intervals that complement the ones passed.

    Examples:
        >>> complement_ranges([(1, 10), (11, 20), (25, 30)], end=100)
        [(21, 24), (31, 100)]
        >>> complement_ranges([(1, 10), (11, 20), (25, 30)])
        [(21, 24)]
        >>> complement_ranges([(0, 2), (3, 17), (18, 20)])
        []
        >>> complement_ranges([(0, 2), (3, 17), (18, 20)], end=100)
        [(21, 100)]
    """
    intervals = union_ranges(intervals)

    comp_intervals = []

    if intervals[0][0] > 1:
        comp_intervals.append((1, intervals[0][0] - 1))

    for index in range(0, len(intervals) - 1):
        new_start = intervals[index][1] + 1
        new_end = intervals[index + 1][0] - 1
        if new_start < new_end:
            comp_intervals.append(
                (new_start, new_end)
            )

    if end is not None:
        start = intervals[-1][1] + 1
        if end > start:
            comp_intervals.append((start, end))

    return comp_intervals


def ranges_length(ranges):
    """
    .. versionadded:: 0.1.12

    Given an iterable where each element is a range, a tuple whose elements
    are numbers with the first being less than or equal to the second, the
    function sums the lengths of all ranges.

    .. warning::

        it's supposed to be used on intervals that were first passed to
        functions like :func:`union_ranges`. If values overlap, there the sum
        will be wrong

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

    for index in range(0, len(data), step):
        yield func(data[index:index+window])


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
            filename=func.__code__.co_filename,
            lineno=func.__code__.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func
