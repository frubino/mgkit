"""
Old GFF filtering module
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


def choose_annotation(ann1, ann2, threshold, only_same_gene=False,
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
                       choose_func=choose_annotation,
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
