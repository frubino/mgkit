"""
Module dealing with BAM/SAM files
"""

import numpy
import logging
import itertools
import pandas

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
        end=norm_end,
        truncate=True
    )

    coverage = pandas.Series(
        dict(
            (pileup_proxy.pos, pileup_proxy.n)
            for pileup_proxy in iterator
        ),
        dtype=numpy.int,
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

    tot_coverage = {}

    for bam_file, sample in itertools.izip(bam_files, samples):
        LOG.info("Adding coverage for sample %s", sample)
        for index, annotation in enumerate(annotations):
            sample_coverage = get_region_coverage(
                bam_file,
                annotation.seq_id,
                annotation.start,
                annotation.end
            )
            #adds the results of the coverage to the total coverage
            #uses the add method of a Series to make sure that possible
            #nan values are filled with 0 before the sum
            try:
                tot_coverage[annotation.uid] = sample_coverage.add(
                    tot_coverage[annotation.uid],
                    fill_value=0
                )
            except KeyError:
                tot_coverage[annotation.uid] = sample_coverage

            cov_mean = sample_coverage.mean()

            annotation.set_attr(
                sample + attr_suff,
                0 if numpy.isnan(cov_mean) else int(cov_mean)
            )
            if (index + 1) % 1000 == 0:
                LOG.debug(
                    'Added coverage for (%d/%d) annotations',
                    index + 1,
                    len(annotations)
                )
    for annotation in annotations:
        cov_mean = tot_coverage[annotation.uid].mean()
        annotation.set_attr(
            'cov',
            0 if numpy.isnan(cov_mean) else int(cov_mean)
        )
