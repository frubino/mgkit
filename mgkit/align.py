"""
Module dealing with BAM/SAM files

.. todo::

    * make a script to add coverage info to gff file
"""

import numpy
import logging
import pandas
from .utils.common import between
import itertools

LOG = logging.getLogger(__name__)


def get_region_coverage(bam_file, seq_id, feat_from, feat_to):
    """
    Return coverage for an annotation.

    .. note::

        feat_from and feat_to are 1-based indexes

    :param Samfile bam_file: instance of :class:`pysam.Samfile`
    :param str seq_id: sequence id
    :param int feat_from: start position of feature
    :param int feat_to: end position of feature

    :return int: coverage array for the annotation
    """
    norm_start, norm_end = feat_from - 1, feat_to - 1
    iterator = bam_file.pileup(
        reference=seq_id,
        start=norm_start,
        end=norm_end
    )

    coverage = pandas.Series(
        dict(
            (pileup_proxy.pos, pileup_proxy.n) for pileup_proxy in iterator
            if between(pileup_proxy.pos, norm_start, norm_end)
        ),
        dtype=numpy.int,
        index=range(norm_start, norm_end + 1)
    )

    return coverage


def add_coverage_info(annotations, bam_files, samples, attr_suff='_cov'):
    """
    Adds coverage information to annotations, using BAM files.

    Thee coverage information is added for each sample as a 'sample_cov' and the
    total coverage as as 'cov' attribute in the annotations.

    .. note::

        The bam_files and sample variables must have the same order

    :param iterable annotations: iterable of annotations
    :param iterable bam_files: iterable of :class:`pysam.Samfile` instances
    :param iterable sample: names of the samples for the BAM files
    """
    LOG.info("Adding coverage info for %d annotations", len(annotations))
    for bam_file, sample in zip(bam_files, samples):
        LOG.info("Sample %s, file %s", sample, bam_file.filename)

    tot_coverage = [
        numpy.zeros(len(annotation))
        for annotation in annotations
    ]

    for bam_file, sample in itertools.izip(bam_files, samples):
        LOG.info("Adding coverage for sample %s", sample)
        for index, annotation in enumerate(annotations):
            sample_coverage = get_region_coverage(
                bam_file,
                annotation.seq_id,
                annotation.feat_from,
                annotation.feat_to
            )
            #adds the results of the coverage to the total coverage
            tot_coverage[index] += sample_coverage
            cov_mean = sample_coverage.mean()
            setattr(
                annotation.attributes,
                sample + attr_suff,
                0 if numpy.isnan(cov_mean) else int(cov_mean)
            )
            if (index + 1) % 1000 == 0:
                LOG.debug(
                    'Added coverage for (%d/%d) annotations',
                    index + 1,
                    len(annotations)
                )
    for index, annotation in enumerate(annotations):
        cov_mean = tot_coverage[index].mean()
        setattr(
            annotation.attributes,
            'cov',
            0 if numpy.isnan(cov_mean) else int(cov_mean)
        )
