"""
GFF filtering
"""
from __future__ import division
import operator
import itertools
from ..utils.common import average_length, between

#------------------------------------------------
#old filter code, most of it should be deprecated


def choose_by_score(ann1, ann2):
    "the winner is the first element of the list"
    return sorted((ann1, ann2), key=operator.attrgetter('score'))


def _choose_annotation(ann1, ann2, threshold, only_same_gene=False,
                       only_same_strand=True, score_func=choose_by_score):
    "Choos one annotation, based on the score"
    #if we don't want to exclude overlaps from different strands
    if only_same_strand:
        if ann1.attributes.frame[0] != ann2.attributes.frame[0]:
            return None, None
    #if we want to filter only genes that have the same id
    if only_same_gene:
        if ann1.attributes.ko != ann2.attributes.ko:
            return None, None

    a1s = ann1.feat_from
    a2s = ann2.feat_from
    a1e = ann1.feat_to
    a2e = ann2.feat_to
    a1s_in_a2 = between(a1s, a2s, a2e)
    a1e_in_a2 = between(a1e, a2s, a2e)
    a2s_in_a1 = between(a2s, a1s, a1e)
    a2e_in_a1 = between(a2e, a1s, a1e)
    a1_a2_avg = average_length(a1s, a1e, a2s, a2e)

    #A1 included in A2 o A2 included in A1
    if (a1s_in_a2 & a1e_in_a2) | (a2s_in_a1 & a2e_in_a1):
        return score_func(ann1, ann2)
    #A1 and A2 overlap
    elif a1e_in_a2 & (a1s < a2s):
        overlap = a1e - a2s + 1
        if overlap / a1_a2_avg > threshold:
            return score_func(ann1, ann2)
    elif a2e_in_a1 & (a2s < a1s):
        overlap = a2e - a1s + 1
        if overlap / a1_a2_avg > threshold:
            return score_func(ann1, ann2)

    #In case no condition was met
    return None, None


def filter_overlapping(annotations, threshold, same_gene, both_strands,
                       choose_func=_choose_annotation,
                       score_func=choose_by_score):
    """
    Filters annotations checking how much they overlap, in which case the ones
    with better score are used.

    :param iterable annotations: iterable of GFF annotations
    :param float threshold: maximum overlap allowed
    :param bool same_gene: if True only filter genes that have the same id
    :param bool both_strands: if True checks annotations on both strands
    :param func choose_func: function used to choose the annotations. Defaults
        to :func:`choose_annotation`
    :param func score_func: function used to decide which annotation to keep
        if two overlap
    :return set: annotations that passed the filtering
    """
    if len(annotations) == 1:
        return annotations

    winners = set()
    losers = set()
    for annotation1, annotation2 in itertools.combinations(annotations, 2):
        winner, loser = choose_func(
            annotation1,
            annotation2,
            threshold,
            same_gene,
            both_strands,
            score_func=score_func
        )
        if winner is None:
            winners.add(annotation1)
            winners.add(annotation2)
        else:
            winners.add(winner)
            losers.add(loser)

    return winners - losers  # [x for x in set(winners) if x not in losers]


#Filters return True if annotation passes a filter and False otherwise
def filter_by_score(threshold, annotation):
    "Filter based on the score of the annotation"
    return True if annotation.score <= threshold else False


def filter_by_bit_score(threshold, annotation):
    "Filter based on the bit score of the annotation"
    return True if float(
        annotation.attributes.bit_score
    ) >= threshold else False


def filter_by_taxon(taxa, annotation):
    "Filter based on the taxon name of the annotation"
    return True if annotation.attributes.taxon.lower() in taxa else False


def filter_by_seq_id(seq_ids, annotation):
    "Filter based on the sequence containing the annotation"
    return True if annotation.seq_id.lower() in seq_ids else False


def filter_by_ko_id(ko_ids, annotation):
    "Filter based on the KO id of the annotation"
    return True if annotation.attributes.ko.lower() in ko_ids else False


def filter_by_ko_idx(ko_ids, annotation):
    "Filter based on the KO id (indexed) of the annotation"
    return True if annotation.attributes.ko_idx.lower() in ko_ids else False


def filter_by_reviewed(annotation):
    "Filter based on the reviewed attribute of the annotation"
    return annotation.attributes.reviewed


def filter_by_strand(strand, annotation):
    "Filter based on the strand of the annotation"
    return True if annotation.strand == strand else False


def filter_by_description(description, annotation):
    "Filter based on the description of the annotation"
    if description in annotation.attributes.description.lower():
        return True
    else:
        return False


def filter_by_hit_length(ko_len, perc, annotation):
    "Filter based on the its length and the profile length of the annotation"
    try:
        avg_len = ko_len[annotation.attributes.name]
    except KeyError:
        avg_len = ko_len[
            (
                annotation.attributes.ko,
                annotation.get_taxon_id(),
                annotation.reviewed
            )
        ]
    hit = (
        int(annotation.attributes.aa_to) -
        int(annotation.attributes.aa_from) + 1
    ) / avg_len
    return hit >= perc


#---------------------------------------


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
