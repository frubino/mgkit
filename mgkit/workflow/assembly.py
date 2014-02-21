"""
Workflow associated with assembly statistics and evaluation
"""

import logging
from ..io import gff
from .. import taxon
from ..utils.sequence import calc_n50
import numpy

LOG = logging.getLogger(__name__)


def rank_annotations_by_attr(annotations, attr='taxon'):
    """
    For all annotations in the list (usually all annotations for a contig),
    counts how many time a set attribute 'attr' appears. The resulting
    dictionary is then sorted by the number of counts and the one with the
    highest count is ranked by how much it represent the total number of
    counts.

    The rank is an integer number between 0 and 10.

    :param iterable annotations: list of :class:`gff.GFFKegg` instances
    :param str attr: the attribute for which the annotations are counted

    :return tuple: the attr with the most counts and its rank
    """
    taxon_counts = {}

    for annotation in annotations:
        taxon_name = getattr(annotation.attributes, attr)
        try:
            taxon_counts[taxon_name] += 1
        except KeyError:
            taxon_counts[taxon_name] = 1

    taxon_counts = sorted(
        (count, taxon_name) for taxon_name, count in taxon_counts.items()
    )
    total_counts = float(sum(x[0] for x in taxon_counts))

    rank = int((taxon_counts[-1][0] / total_counts) * 10)

    return taxon_counts[-1][1], rank


def assign_contigs_to_taxa(annotations, root_map=None,
                           black_list=None):
    """
    Groups annotations by contig (seq_id) and counts how many contigs a taxon,
    or its root if root_map is supplied, have been assigned to.

    The actual form of the dictionary like this:

    .. digraph:: results
        :alt: taxon->ranks->count

        "taxon 1" -> "rank 0" -> count;
        "taxon 1" -> "rank 1" -> count;
        "taxon 1" -> "rank 4" -> count;
        "taxon 1" -> "rank 7" -> count;
        "taxon 1" -> "rank 9" -> count;
        "taxon 1" -> "rank 10" -> count;

    .. note::

        the number of ranks for a taxon is not pretedermined, but depends on
        the values returned by :func:`rank_annotations_by_attr`.

    :param iterable annotations: list of :class:`gff.GFFKegg` instances
    :param dict root_map: dictionary taxon->root

    :return dict: dictionary
    """

    LOG.info("Assigning %d annotations", len(annotations))
    LOG.info(
        "Blacklisted taxa: %s",
        'none' if root_map is None else ', '.join(black_list)
    )
    if root_map is not None:
        LOG.info("Using only root taxa")

    contigs = gff.group_by_contig(annotations)

    assignments = {}
    for contig, ann_list in contigs.iteritems():
        taxon_name, rank = rank_annotations_by_attr(ann_list, attr='taxon')

        if black_list is not None:
            if taxon_name in black_list:
                continue

        if root_map is not None:
            taxon_name = taxon.get_taxon_root(taxon_name, root_map)

        if taxon_name not in assignments:
            assignments[taxon_name] = {}
        value = 1
        try:
            assignments[taxon_name][rank] += value
        except KeyError:
            assignments[taxon_name][rank] = value

    return assignments


def filter_contig_assignments(contig_assign, threshold=5, min_counts=1):
    """
    Filter contigs assignments using a threshold for the rank: all rank counts
    belonging to a taxon which are greater than or equal to threshold will be
    summed up.

    :param dict contig_assign: dictionary returned by
        :func:`assign_contigs_to_taxa`
    :param int threshold: the minimum rank for which the counts are summed up

    :return dict: dictionary in the form taxon_name->count
    """

    LOG.info(
        "Filtering contig assignments with a minimum rank of %d and a " +
        "minimum count of %d", threshold, min_counts
    )

    new_counts = {}
    for taxon, rank_dict in contig_assign.items():

        rank_sum = sum(
            count for rank, count in rank_dict.items()
            if rank >= threshold
        )

        if rank_sum < min_counts:
            continue

        new_counts[taxon] = rank_sum

    return new_counts


def basic_stats(array, sep):
    """
    Returns formatted basic statistics for contig lengths
    """
    amedian = numpy.median(array)
    amean = numpy.mean(array)
    amax = numpy.max(array)
    amin = numpy.min(array)
    return sep.join(str(value) for value in (amedian, amean, amax, amin))


def write_fasta_summary(file_handle, seq_lengths, seq_lengths_filt, sep='\t'):
    """
    Write summary file for assembly

    :param file_handle: file handle for output
    :param array seq_lengths: array for sequence lengths
    :param array seq_lengths_filt: array for sequence lengths of annotated
        contigs
    :param str sep: string used as column separator
    """
    n50 = calc_n50(seq_lengths)
    file_handle.write("Assembly N50:{0}{1}\n\n".format(sep, n50))
    header = "{0}\n".format(sep.join(("", "median", "mean", "max", "min")))
    file_handle.write(header)
    file_handle.write("assembly{0}{1}\n".format(sep,
                      basic_stats(seq_lengths, sep)))
    file_handle.write("annotated contigs{0}{1}\n".format(sep,
                      basic_stats(seq_lengths_filt, sep)))
