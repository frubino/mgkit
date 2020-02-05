"""
Module dealing with BAM/SAM files
"""
from future.utils import viewitems
import logging
import weakref
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
    .. versionchanged:: 0.4.2
        the function returns **lists** instead of numpy arrays for speed (at
        least in my tests it seems ~4x increase)

    .. versionchanged:: 0.4.0
        now returns 3 array, instead of 2. Also added *seq_ids* to skip lines

    .. versionchanged:: 0.3.4
        *num_seqs* can be None to avoid a log message

    .. versionadded:: 0.3.0

    Reads a samtools *depth* file, returning a generator that yields the
    array of each base coverage on a per-sequence base.

    .. note::

        There's no need anymore to use `samtools depth -aa`, because the
        function returns the position array and this can be used to create a
        Pandas SparseArray which can be reindexed to include missing positions
        (with values of 0)

        **Valid for version < 0.4.0**:

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
        is the list with the positions, the third element is the list with the
        coverages
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
        # From Python3 the default is Universal newlines, and it's not expected
        # to have more than '\n' at the end of the line - increases speed
        # slightly
        name, pos, cov = line[:-1].split('\t')
        if (seq_ids is not None) and (name not in seq_ids):
            continue
        # only converts if sequence is to be used
        pos = int(pos)
        cov = int(cov)

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
                yield curr_key, curr_pos, curr_cov
                curr_key = name
                curr_cov = [cov]
                curr_cov = [pos]
    else:
        yield curr_key, curr_pos, curr_cov

    LOG.info('Read a total of %d sequence coverage', line_no + 1)


class SamtoolsDepth(object):
    """
    .. versionchanged:: 0.4.2
        several optimisations and changes to support a scanning approach,
        instead of lookup table. No exception is raised when a sequence is not
        found in the file, instead assuming that the coverage is 0

    .. versionchanged:: 0.4.0
        uses pandas.SparseArray now. It should use less memory, but needs
        pandas version > 0.24

    .. versionadded:: 0.3.0

    This a class that helps use the results from :func:`read_samtools_depth`.

    There are 2 modes of operations:

        1) Request a region coverage via :meth:`SamtoolsDepth.region_coverage`

        2) Advance the `samtools depth` file until a sequence coverage is read

    The method 1) was the default in MGKit < 0.4.2, when the user would request
    a region coverage and this class would advance the reading until the
    requested region is found. While scanning, all the sequences encountered
    are kept as pandas Series with SparseArray declared, reindexed from 0 to
    the *max_size_dict* values if provided, or *max_size*. The advantage is that
    the it's easier to use, but with bigger datasets will keep high memory
    usage. This is the case for case 2) which involves calling directly
    :meth:`SamtoolsDepth.advance_file` returning the sequence ID read. Then the
    region coverage can be requested the same way as before, but won't involve
    scanning the file.

    The use of a SparseArray in a pandas Series allows the use of `samtools depth` files
    that weren't produced with `samtools depth -aa`.

    .. note::

        Starting with MGKit 0.4.2, the internal dictionary to keep the SparseArray(s)
        is a :class:`weakref.WeakValueDictionary`, which should improve the release
        of memory. However, the amount of memory used is still fairly high, especially
        with the increasing number of sequences in a GFF/Depth file. It is recoommended
        to use max_size_dict to 1) only creates arrays for sequences needed and 2) use
        arrays of smaller size

    """
    file_handle = None
    data = None
    max_size = None
    max_size_dict = None
    not_found = None
    closed = False
    density = None

    def __init__(self, file_handle, num_seqs=10**4, max_size=10**6,
                 max_size_dict=None, calc_density=False, dtype='uint32'):
        """
        .. versionchanged:: 0.4.2
            added *raise_error* to control sequences not found in depth files,
            *calc_density* to debug the density of the SparseArray used and
            also dtype to control the dtype to use in the SparseArray

        .. versionchanged:: 0.4.0
            added *max_size* and max_size_dict

        Arguments:
            file_handle (file): the file handle to pass to
                :func:`read_samtools_depth`
            num_seqs (int): number of sequence that fires a log message
            max_size (int): max size to use in the SparseArray
            max_size_dict (dict): dictionary with max size for each seq_id
                requested. If None, *max_size* is used
            calc_density (bool): If True, the density of the SparseArray(s)
                used are kept in :attr:`SamtoolsDepth.density`
            dtype (str): dtype to pass to the SparseArray

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
        self.not_found = set()
        self.closed = False
        self.data = weakref.WeakValueDictionary()
        self.dtype = 'Sparse[{}]'.format(dtype)
        if calc_density:
            self.density = []

    def advance_file(self):
        """
        .. versionadded:: 0.4.2

        Requests the scan of the file and returns the sequence found, the
        coverage can then be retrivied with :meth:`SamtoolsDepth.region_coverage`

        Returns:
            (str, None): the return value is None if the file is fully read,
            otherwise the the sequnece ID found is returned
        """
        # continue scanning the file
        try:
            key, pos_array, cov_array = next(self.file_handle)
        except StopIteration:
            self.closed = True
            return None
        # Init a SparseArray (by passing the correct dtype)
        value = pandas.Series(
            # brings start position to 0
            dict(zip((pos - 1 for pos in pos_array), cov_array)),
            dtype=self.dtype
        # reindex to the correct size to allow
        ).reindex(
            range(self.max_size_dict.get(key, self.max_size)),
            fill_value=0
        )
        self.data[key] = value
        if self.density is not None:
            self.density.append(value.sparse.density)
        return key

    def region_coverage(self, seq_id, start, end):
        """
        .. versionchanged:: 0.4.2
            now using :meth:`SamtoolsDepth.advance_file` to scan the file

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
        if seq_id in self.data:
            # data already in dictionary
            return self.data[seq_id][start-1:end].mean()
        elif seq_id in self.not_found:
            # already asked for and not found in depth file
            return 0.
            # depth file is closed and we assume
            # so no coverage is available
        elif self.closed:
            self.not_found.add(seq_id)
            return 0.
        else:
            while True:
                key = self.advance_file()
                if key is None:
                    self.not_found.add(seq_id)
                    self.closed = True
                    return 0.
                elif key == seq_id:
                    return self.data[seq_id][start-1:end].mean()

    def drop_sequence(self, seq_id):
        """
        .. versionadded:: 0.4.2

        Remove the sequence passed from the internal dictionary
        """
        if seq_id in self.data:
            del self.data[seq_id]
