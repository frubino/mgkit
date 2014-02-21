"""
Module that contains mapping operations on pandas data structures
"""

import pandas


def group_dataframe_by_mapping(dataframe, mapping, root_taxon, name_dict=None):
    """
    Return a :class:`pandas.DataFrame` filtered by mapping and root taxon, the
    values for each column is averaged over all genes mapping to a category.

    :param DataFrame dataframe: DataFrame with multindex gene-root
    :param dict mapping: dictionary of category->genes
    :param str root_taxon: root taxon to group genes
    :param dict name_dict: dictionary of category->name

    :return DataFrame: DataFrame filtered
    """
    data_dict = {}
    for category, id_list in mapping.iteritems():
        category_name = category
        if name_dict is not None:
            category_name = name_dict[category_name]
        data_dict[category_name] = dataframe.loc[id_list].xs(
            root_taxon, level='root'
        ).mean()

    return pandas.DataFrame.from_dict(data_dict, orient='index')


def calc_coefficient_of_variation(dataframe):
    """
    Calculate coefficient of variation for a DataFrame. Uses formula from
    `Wikipedia <http://en.wikipedia.org/wiki/Coefficient_of_variation>`_

    The formula used is :math:`\\left (1 + \\frac {1}{4n} \\right ) * c_{v}`
    where :math:`c_{v} = \\frac {s}{\\bar{x}}`
    """

    dataframe_cv = dataframe.std(axis=1) / dataframe.mean(axis=1)
    coeff = 1 + (1.0 / (4 * dataframe.count(axis=1)))

    return coeff * dataframe_cv


def make_stat_table(dataframes, roots):
    """
    Produces a :class:`pandas.DataFrame` that summarise the supplied DataFrames.
    The stats include mean, stdev and coefficient of variation for each root
    taxon.

    :param iterable dataframes: iterable of DataFrame instances
    :param iterable roots: list of root taxa to which each table belongs

    :return DataFrame: returns a DataFrame instance
    """

    index = []

    data = {}

    for dataframe, root in zip(dataframes, roots):
        index.append((root, 'mean'))
        data[(root, 'mean')] = dataframe.mean(axis=1)
        index.append((root, 'stdev'))
        data[(root, 'stdev')] = dataframe.std(axis=1)
        index.append((root, 'c. var'))
        data[(root, 'c. var')] = calc_coefficient_of_variation(dataframe)

    index = pandas.MultiIndex.from_tuples(sorted(index),
                                          names=('root', 'value'))

    return pandas.DataFrame(data, columns=index)


def concatenate_and_rename_tables(dataframes, roots):
    """
    Concatenates a list of :class:`pandas.DataFrame` instances and renames the
    columns prepending a string to each column in each table from a list of
    prefixes.

    :param iterable dataframes: iterable of DataFrame instances
    :param iterable roots: list of prefixes to append to the column names of
        each DataFrame

    :return DataFrame: returns a DataFrame instance

    .. todo::

        * move to pandas_utils?
    """
    renamed = []

    for dataframe, root in zip(dataframes, roots):
        index = dict((column, root + column) for column in dataframe.columns)
        renamed.append(dataframe.rename(columns=index))

    return pandas.concat(renamed, axis=1)
