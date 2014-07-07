"""
GFF filtering
"""
from __future__ import division
import itertools


def choose_annotation(ann1, ann2, overlap=100, choose_func=None):
    """
    .. versionadded:: 0.1.12

    Given two :class:`mgkit.io.gff.Annotation`, if one of of the two annotations
    either is contained in the other or they overlap for at least a *overlap*
    number of bases, *choose_func* will be applied to both. The result of
    *choose_func* is the the annotation to be discarderd. It returns *None* if
    the annotations should be both kept.

    No checks are made to ensure that the two annotations are on the same
        sequence and strand, as the *intersect* method of
        :class:`mgkit.io.gff.Annotation` takes care of them.

    Arguments:
        ann1: instance of :class:`mgkit.io.gff.Annotation`
        ann2: instance of :class:`mgkit.io.gff.Annotation`
        overlap (int, float): number of bases overlap that trigger the filtering
        choose_func (None, func): function that accepts *ann1* and *ann2* and
            return the one to be discarded or None if both are accepted

    Returns:
        (None, Annotation): returns either the :class:`mgkit.io.gff.Annotation`
        to be discarded or None, which is the result of *choose_func*

    .. note::

        If *choose_func* is *None*, the default function is used::

            lambda a1, a2: min(a1, a2, key=lambda el: (el.dbq, el.bitscore, len(el)))

        In order of importance the db quality, the bitscore and the length. The
        annotation with the lowest tuple value is the one to discard.

    """

    if choose_func is None:
        choose_func = lambda a1, a2: min(a1, a2, key=lambda el: (el.dbq, el.bitscore, len(el)))

    intersect = ann1.intersect(ann2)

    if intersect is not None:
        #if the intersection is the same size of one of the annotations size,
        #it means that one contain the other
        if (len(intersect) == len(ann1)) or (len(intersect) == len(ann2)):
            return choose_func(ann1, ann2)
        else:
            #if the overlap is longer than the threshold
            if len(intersect) > overlap:
                return choose_func(ann1, ann2)

    return None


def filter_annotations(annotations, choose_func=None, sort_func=None, reverse=True):
    """
    .. versionadded:: 0.1.12

    Filter an iterable of :class:`mgkit.io.gff.Annotation` instances sorted
    using *sort_func* as key in *sorted* and if the order is to be *reverse*;
    it then applies *choose_func* on all possible pair combinations, using
    itertools.combinations.

    By default *choose_func* is :func:`choose_annotation` with the default
    values, the list of annotation is sorted by bitscore, from the highest to
    the lowest value.

    Arguments:
        annotations (iterable): iterable of :class:`mgkit.io.gff.Annotation`
            instances
        choose_func (func, None): function used to select the *losing*
            annotation; if `None`, it will be :func:`choose_annotation` with
            default values
        sort_func (func, None): by default the sorting key is the bitscore of
            the annotations
        reverse (bool): passed to `sorted`, by default is reversed

    Returns:
        set: a set with the annotations that pass the filtering
    """

    if choose_func is None:
        choose_func = choose_annotation

    if sort_func is None:
        sort_func = lambda x: x.bitscore

    annotations = sorted(annotations, key=sort_func, reverse=reverse)

    to_remove = set()

    for ann1, ann2 in itertools.combinations(annotations, 2):

        if ann1 in to_remove or ann2 in to_remove:
            continue

        to_remove.add(choose_func(ann1, ann2))

    return set(annotations) - to_remove


def filter_base(annotation, attr=None, value=None):
    """
    Checks if an annotation attribute is equal to the requested value

    Arguments:
        annotation: :class:`mgkit.io.gff.Annotation` instance
        attr (str): attribute of the annotation
        value: the value that the attribute should be equal to

    Returns:
        bool: True if the supplied value is equal to the attribute ot False
        otherwise
    """
    return getattr(annotation, attr, None) == value


def filter_len(annotation, value=None, greater=True):
    """
    Checks if an annotation length is longer, equal of shorter than the
    requested value

    Arguments:
        annotation: :class:`mgkit.io.gff.Annotation` instance
        value (int): the length to which the attribute should be compared to
        greater (bool): if True the annotation length must be equal or greater
            than and if False equal of lower than

    Returns:
        bool: True if the test passes
    """
    if greater:
        return len(annotation) >= value
    else:
        return len(annotation) <= value


def filter_base_num(annotation, attr=None, value=None, greater=True):
    """
    Checks if an annotation attribute is greater, equal of lower than the
    requested value

    Arguments:
        annotation: :class:`mgkit.io.gff.Annotation` instance
        attr (str): attribute of the annotation
        value (int): the value to which the attribute should be compared to
        greater (bool): if True the attribute value must be equal or greater
            than and if False equal of lower than

    Returns:
        bool: True if the test passes
    """
    annotation_value = getattr(annotation, attr, None)
    if annotation_value is None:
        return False
    annotation_value = float(annotation_value)
    if greater:
        return annotation_value >= value
    else:
        return annotation_value <= value


def filter_attr_num(annotation, attr=None, value=None, greater=True):
    """
    Checks if an annotation *attr* dictionary contains a key shose value is
    greater, equal of lower than the requested value

    Arguments:
        annotation: :class:`mgkit.io.gff.Annotation` instance
        attr (str): key in the :attr:`mgkit.io.gff.Annotation.attr`
            dictionary
        value (int): the value to which we need to compare
        greater (bool): if True the value must be equal or greater than and if
            False equal of lower than

    Returns:
        bool: True if the test passes
    """
    annotation_value = annotation.attr.get(attr, None)
    if annotation_value is None:
        return False
    annotation_value = float(annotation_value)
    if greater:
        return annotation_value >= value
    else:
        return annotation_value <= value


def filter_attr_str(annotation, attr=None, value=None, equal=True):
    """
    Checks if an annotation *attr* dictionary contains a key shose value is
    equal to, or contains the requested value

    Arguments:
        annotation: :class:`mgkit.io.gff.Annotation` instance
        attr (str): key in the :attr:`mgkit.io.gff.Annotation.attr`
            dictionary
        value (int): the value to which we need to compare
        equal (bool): if True the value must be equal and if False equal value
            must be contained

    Returns:
        bool: True if the test passes
    """
    annotation_value = annotation.attr.get(attr, None)
    if annotation_value is None:
        return False
    if equal:
        return annotation_value == value
    else:
        return value in annotation_value
