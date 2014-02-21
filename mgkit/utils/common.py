"""
Utility functions
"""

import warnings


def avg_len(a1s, a1e, a2s, a2e):
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


def deprecated(func):
    '''
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.

    .. todo::

        needs to check the correct license or rewrite

    '''
    def new_func(*args, **kwargs):
        "Wrapper for the function"
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func
