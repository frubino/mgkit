"""
Utilities to map genes
"""

import pandas
import itertools


def count_genes_in_mapping(gene_lists, labels, mapping, normalise=False):
    """
    Maps lists of ids to a mapping dictionary, returning a
    :class:`pandas.DataFrame` in which the rows are the labels provided and
    the columns the categories to which the ids map. Each element of the matrix
    label-category is the sum of all ids in the relative gene list that maps to
    the specific category.

    :param iterable gene_lists: an iterable in which each element is a iterable
        of ids that can be mapped to mapping
    :param iterable labels: an iterable of strings that defines the labels to
        be used in the resulting rows in the :class:`pandas.DataFrame`; must
        have the same length as gene_lists
    :param dict mapping: a dictionary in the form:
        gene_id->[cat1, cat2, .., catN]
    :param bool normalise: if True the counts are normalised over the total for
        each row.

    :return: a :class:`pandas.DataFrame` instance
    """
    categories = set(itertools.chain(*mapping.values()))

    matrix = pandas.DataFrame(index=labels, columns=categories)
    matrix.fillna(0, inplace=True)

    for label, gene_list in zip(labels, gene_lists):
        for gene_id in gene_list:
            try:
                mapped_ids = mapping[gene_id]
            except KeyError:
                continue

            for mapped_id in mapped_ids:
                matrix.ix[label, mapped_id] += 1

    matrix = matrix.ix[:, matrix.sum() > 0]

    if normalise:
        matrix = matrix.div(matrix.sum(axis=1), axis=0)

    return matrix


def group_annotation_by_mapping(annotations, mapping, attr='ko'):
    """
    Group annotations by mapping dictionary

    :param iterable annotations: iterable of :class:`gff.GFFKeg` instances
    :param dict mapping: dictionary with mappings for the attribute requested
    :param str attr: attribute of the annotation to be used as key in mapping

    :return dict: dictionary category->annotations
    """
    grouped = dict(
        (categ, []) for categ in set(itertools.chain(*mapping.values()))
    )
    for annotation in annotations:
        try:
            categories = mapping[getattr(annotation.attributes, attr)]
        except KeyError:
            #not included in the categories
            continue
        for categ in categories:
            grouped[categ].append(annotation)

    return grouped
