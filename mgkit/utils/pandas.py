"""
Functions related to handling :mod:`pandas` data structures
"""


def filter_dataframe_by_count(dataframe, min_count, no_empty_rows=True):
    """
    Filter the number of column of a DataFrame based on a minimum amount of
    non-NA counts and removes all rows that have no non-NA values if
    no_empty_rows is True (default)

    :param DataFrame dataframe: :class:`pandas.DataFrame` instance
    :param int min_count: minimum number of non-NA counts for columns
    :param bool no_empty_rows: if True, removes all empty rows

    :return: :class:`pandas.DataFrame` instance
    """

    cols = dataframe.ix[:, dataframe.count(axis=0) > min_count].columns
    dataframe = dataframe.ix[:, cols]
    if no_empty_rows:
        rows = dataframe[dataframe.count(axis=1) > 0].index
    else:
        rows = dataframe.index

    return dataframe.ix[rows]
