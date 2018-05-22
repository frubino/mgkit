"""
Module dealing with BAM/SAM files
"""
from future.utils import viewitems
import logging
import itertools
try:
    # In Python2
    from itertools import izip as zip
except ImportError:
    pass
from builtins import object
import pandas
import numpy
import pysam
import progressbar
from mgkit.io.utils import open_file

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


def covered_annotation_bp(files, annotations, min_cov=1, progress=False):
    """
    .. versionadded:: 0.1.14

    Returns the number of base pairs covered of annotations over multiple
    samples.

    Arguments:
        files (iterable): an iterable that returns the alignment file names
        annotations (iterable): an iterable that returns annotations
        min_cov (int): minumum coverage for a base to counted
        progress (bool): if *True*, a progress bar is used

    Returns:
        dict: a dictionary whose keys are the uid and the values the number of
        bases that are covered by reads among all samples

    """
    annotations = [
        (annotation.uid, annotation.seq_id, annotation.start, annotation.end)
        for annotation in annotations
    ]

    covered = {}

    for file_name in files:
        # pysam 0.8.1 changed Samfile to AlignmentFile
        alg_file = pysam.AlignmentFile(file_name, 'rb')

        if progress:
            bar = progressbar.ProgressBar(max_value=len(annotations))
            annotations = bar(annotations)

        for uid, seq_id, start, end in annotations:
            cov = get_region_coverage(alg_file, seq_id, start, end)

            try:
                covered[uid] = cov.add(covered[uid], fill_value=0)
            except KeyError:
                covered[uid] = cov

    return dict(
        (uid, len(cov[cov >= min_cov])) for uid, cov in viewitems(covered)
    )


def add_coverage_info(annotations, bam_files, samples, attr_suff='_cov'):
    """
    .. versionchanged:: 0.3.4
        the coverage now is returned as floats instead of int

    Adds coverage information to annotations, using BAM files.

    The coverage information is added for each sample as a 'sample_cov' and the
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

    for bam_file, sample in zip(bam_files, samples):
        LOG.info("Adding coverage for sample %s", sample)
        for index, annotation in enumerate(annotations):
            sample_coverage = get_region_coverage(
                bam_file,
                annotation.seq_id,
                annotation.start,
                annotation.end
            )
            # adds the results of the coverage to the total coverage
            # uses the add method of a Series to make sure that possible
            # nan values are filled with 0 before the sum
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
            0 if numpy.isnan(cov_mean) else cov_mean
        )


def read_samtools_depth(file_handle, num_seqs=10000):
    """
    ..versionchanged:: 0.3.4
        *num_seqs* can be None to avoid a log message

    .. versionadded:: 0.3.0

    Reads a samtools *depth* file, returning a generator that yields the
    array of each base coverage on a per-sequence base.

    .. note::

        The information on position is not used, to use numpy and save memory.
        samtools *depth* should be called with the `-aa` option::

             `samtools depth -aa bamfile`

        This options will output both base position with 0 coverage and
        sequneces with no aligned reads

    Arguments:
        file_handle (file): file handle of the coverage file
        num_seqs (int or None): number of sequence that fires a log message. If
            None, no message is triggered

    Yields:
        tuple: the first element is the sequence identifier and the second one
        is the *numpy* array with the positions
    """
    curr_key = ''
    curr_cov = []

    file_handle = open_file(file_handle, 'rb')

    LOG.info(
        'Reading coverage from file (%s)',
        getattr(file_handle, 'name', repr(file_handle))
    )
    line_no = 0
    for line in file_handle:
        line = line.decode('ascii')
        name, pos, cov = line.strip().split('\t')
        cov = int(cov)
        if curr_key == name:
                curr_cov.append(cov)
        else:
            if curr_key == '':
                curr_cov.append(cov)
                curr_key = name
            else:
                line_no += 1
                if (num_seqs is not None) and (line_no % num_seqs == 0):
                    LOG.info('Read %d sequence coverage', line_no)
                yield curr_key, numpy.array(curr_cov)
                curr_key = name
                curr_cov = [cov]
    else:
        yield curr_key, numpy.array(curr_cov)

    LOG.info('Read a total of %d sequence coverage', line_no + 1)


class SamtoolsDepth(object):
    """
    .. versionadded:: 0.3.0

    A class used to cache the results of :func:`read_samtools_depth`, while
    reading only the necessary data from a`samtools depth -aa` file.
    """
    file_handle = None
    data = None

    def __init__(self, file_handle, num_seqs=10**4):
        """
        Arguments:
            file_handle (file): the file handle to pass to
                :func:`read_samtools_depth`
            num_seqs (int): number of sequence that fires a log message
        """
        self.file_handle = read_samtools_depth(file_handle, num_seqs=num_seqs)
        self.data = {}

    def region_coverage(self, seq_id, start, end):
        """
        Returns the mean coverage of a region. The *start* and *end* parameters
        are expected to be 1-based coordinates, like the correspondent
        attributes in :class:`mgkit.io.gff.Annotation` or
        :class:`mgkit.io.gff.GenomicRange`.

        If the sequence for which the coverage is requested is not found, the
        *depth* file is read (and cached) until it is found.

        Arguments:
            seq_id (str): sequence for which to return mean coverage
            start (int): start of the region
            end (int): end of the region

        Returns:
            float: mean coverage of the requested region
        """
        try:
            cov = self.data[seq_id][start-1:end]
        except KeyError:
            for key, value in self.file_handle:
                self.data[key] = value
                # if the key is the one requested, the loop is stopped and and
                # the value is kept
                if key == seq_id:
                    cov = value[start-1:end]
                    break

        return cov.mean()
