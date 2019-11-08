"""
Module dealing with BAM/SAM files
"""
from future.utils import viewitems
import logging
try:
    # In Python2
    from itertools import izip as zip
except ImportError:
    pass
from builtins import object
import pandas
import numpy
import pysam
from tqdm import tqdm
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
            annotations = tqdm(annotations)

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


def read_samtools_depth(file_handle, num_seqs=10000, seq_ids=None):
    """
    .. versionchanged:: 0.4.0
        now returns 3 array, instead of 2. Also added *seq_ids* to skip lines

    .. versionchanged:: 0.3.4
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
        seq_ids (dict, set): a hashed container like a dictionary or set with
            the sequences to return

    Yields:
        tuple: the first element is the sequence identifier, the second one
        is the *numpy* array with the positions, the third element is the
        *numpy* array with the coverages
    """
    curr_key = ''
    curr_pos = []
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
        pos = int(pos)
        cov = int(cov)
        if (seq_ids is not None) and (name not in seq_ids):
            continue
        if curr_key == name:
                curr_pos.append(pos)
                curr_cov.append(cov)
        else:
            if curr_key == '':
                curr_cov.append(cov)
                curr_pos.append(pos)
                curr_key = name
            else:
                line_no += 1
                if (num_seqs is not None) and (line_no % num_seqs == 0):
                    LOG.info('Read %d sequence coverage', line_no)
                yield curr_key, numpy.array(curr_pos), numpy.array(curr_cov)
                curr_key = name
                curr_cov = [cov]
                curr_cov = [pos]
    else:
        yield curr_key, numpy.array(curr_pos), numpy.array(curr_cov)

    LOG.info('Read a total of %d sequence coverage', line_no + 1)


class SamtoolsDepth(object):
    """
    .. versionchanged:: 0.4.0
        uses pandas.SparseArray now. It should use less memory, but needs
        pandas version > 0.24

    .. versionadded:: 0.3.0

    A class used to cache the results of :func:`read_samtools_depth`, while
    reading only the necessary data from a`samtools depth -aa` file.
    """
    file_handle = None
    data = None
    max_size = None
    max_size_dict = None

    def __init__(self, file_handle, num_seqs=10**4, max_size=10**6,
                 max_size_dict=None):
        """
        .. versionchanged:: 0.4.0
            added *max_size* and max_size_dict

        Arguments:
            file_handle (file): the file handle to pass to
                :func:`read_samtools_depth`
            num_seqs (int): number of sequence that fires a log message
            max_size (int): max size to use in the SparseArray
            max_size_dict (dict): dictionary with max size for each seq_id
                requested. If None, *max_size* is used
        """
        self.max_size = max_size
        if max_size_dict is None:
            self.max_size_dict = {}
        else:
            self.max_size_dict = max_size_dict

        self.file_handle = read_samtools_depth(
            file_handle,
            num_seqs=num_seqs,
            seq_ids=max_size_dict
        )
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
            for key, pos_array, cov_array in self.file_handle:
                # If the max_size_dict was passed, skips arrays that are not in
                # of interests
                if self.max_size_dict and (key not in self.max_size_dict):
                    continue
                # Init a SparseArray (by passing the correct dtype)
                value = pandas.Series(
                    # brings start position to 0
                    dict(zip(pos_array - 1, cov_array)),
                    dtype="Sparse[int]"
                # reindex to the correct size to allow
                ).reindex(
                    range(self.max_size_dict.get(seq_id, self.max_size)),
                    fill_value=0
                )
                self.data[key] = value
                # if the key is the one requested, the loop is stopped and and
                # the value is kept
                if key == seq_id:
                    cov = value[start-1:end]
                    break
            else:
                raise ValueError(
                    "No coverage information found for {}".format(seq_id)
                )

        return cov.mean()
