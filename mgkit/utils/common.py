"""
Utility functions
"""
from __future__ import division


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

    Given an iterable of where each element is a range, a tuple whose elements
    are numbers with the first being less than or equal to the second, the
    function sums the lengths of all ranges.

    Arguments:
        ranges (iterable): each element is a tuple like *(1, 10)*

    Returns:
        int: sum of all ranges lengths

    """
    return sum(range[1] - range[0] + 1 for range in ranges)
