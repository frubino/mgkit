"""
Module used to filter lists
"""

import numpy

from ..consts import MIN_COV


def aggr_filtered_list(val_list, aggr_func=numpy.mean,
                       filt_func=lambda x: False if x < MIN_COV else True):
    """
    Aggregate a list of values using 'aggr_func' on a list that passed the
    filtering in 'filt_func'.

    'filt_func' is a function that returns True or False for each value in
    val_list. If the return value is True, the element is included in the
    values passed to 'aggr_func'. Internally a list comprehension is used and
    the result passed to 'aggr_func'

    :param iterable val_list: list of values
    :param func aggr_func: function used to aggregate the list values
    :param func filt_func: function the return True or False

    :return: the result of the applied 'aggr_func'
    """
    return aggr_func([val for val in val_list if filt_func(val)])
