"""
This modules define classes and function related to manipulation of GFF/GTF
files.

Needs to be adapted for a more general GFF dialect. It incorporates data from
metagenomic analysis right now. BaseGFF is the class to use for a more general
approach.
"""
from __future__ import print_function
from __future__ import division

import random
import itertools
import logging
import uuid
import urllib
import mgkit.io
from ..utils import sequence as seq_utils
from ..consts import MIN_COV
from ..utils.common import between, union_range, ranges_length
from .. import taxon

LOG = logging.getLogger(__name__)


class AttributeNotFound(Exception):
    """
    Raised if an attribute is not found in a GFF file
    """
    pass


def write_gff(annotations, file_handle, verbose=True):
    """
    .. versionchanged:: 0.1.12
        added *verbose* argument

    Write a GFF to file

    Arguments:
        annotations (iterable): iterable that returns :class:`GFFKegg`
            or :class:`Annotation` instances
        file_handle (str, file): file name or file handle to write to
        verbose (bool): if True, a message is logged
    """

    if isinstance(file_handle, str):
        file_handle = open(file_handle, 'w')

    if verbose:
        LOG.info(
            "Writing annotations to file (%s)",
            getattr(file_handle, 'name', repr(file_handle))
        )

    for annotation in annotations:
        annotation.to_file(file_handle)


class GenomicRange(object):
    """
    Defines a genomic range
    """
    seq_id = 'None'
    "Sequence ID"
    strand = '+'
    "Strand"
    start = None
    "Start of the range, 1-based"
    end = None
    "End of the range 1-based"

    def __init__(self, seq_id='None', start=1, end=1, strand='+'):
        self.seq_id = seq_id
        self.strand = strand
        self.start = start
        self.end = end

    def __eq__(self, other):
        if (self.seq_id != other.seq_id) or (self.strand != other.strand):
            return False
        if (self.start != other.start) or (self.end != other.end):
            return False
        return True

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return "{0}({1}):{2}-{3}".format(
            self.seq_id,
            self.strand,
            self.start,
            self.end
        )

    def __repr__(self):
        return str(self)

    def union(self, other):
        """
        Return the union of two :class:`GenomicRange`
        """
        if (self.seq_id == other.seq_id) and (self.strand == other.strand):
            result = union_range(self.start, self.end, other.start, other.end)
            if result is not None:

                gen_range = GenomicRange()
                gen_range.seq_id = self.seq_id
                gen_range.strand = self.strand
                gen_range.start = result[0]
                gen_range.end = result[1]

                return gen_range

        return None

    def expand_from_list(self, others):
        """
        Expand the :class:`GenomicRange` range instance with a list of
        :class:`GenomicRange`

        Arguments:
            others (iterable): iterable of :class:`GenomicRange`
        """
        new_range = self

        for other in others:
            union = new_range.union(other)
            if union is None:
                continue
            new_range = union

        self.start = new_range.start
        self.end = new_range.end

    def intersect(self, other):
        """
        Return an instance of :class:`GenomicRange` that represent the
        intersection of the current instance and another.
        """
        if (self.seq_id == other.seq_id) and (self.strand == other.strand):

            if between(other.start, self.start, self.end) or \
               between(other.end, self.start, self.end) or \
               between(self.start, other.start, other.end) or \
               between(self.end, other.start, other.end):

                gen_range = GenomicRange()
                gen_range.start = max(self.start, other.start)
                gen_range.end = min(self.end, other.end)

                return gen_range

        return None

    def get_range(self):
        """
        .. versionadded:: 0.1.13

        Returns the start and end position as a tuple
        """
        return (self.start, self.end)


class Annotation(GenomicRange):
    """
    .. versionadded:: 0.1.12

    Alternative implementation for an Annotation. When initialised, If *uid* is
    None, a unique id is added using `uuid.uuid4`.
    """
    source = 'None'
    "Annotation source"
    feat_type = 'None'
    "Annotation type (e.g. CDS, gene, exon, etc.)"
    score = 0.0
    "Score associated to the annotation"
    phase = 0
    "Annotation phase, (0, 1, 2)"
    attr = None
    "Dictionary with the key value pairs in the last column of a GFF/GTF"

    def __init__(self, seq_id='None', start=1, end=1, strand='+', source='None', feat_type='None', score=0.0, phase=0, uid=None, **kwd):
        super(Annotation, self).__init__(
            seq_id=seq_id,
            start=start,
            end=end,
            strand=strand
        )

        self.source = source
        self.feat_type = feat_type
        self.score = score
        self.phase = phase

        self.attr = kwd

        if uid is None:
            self.uid = str(uuid.uuid4())
        else:
            self.uid = uid

    def get_ec(self, level=4):
        """
        .. versionadded:: 0.1.13

        Returns the EC values associated with the annotation, cutting them at
        the desired level.

        Arguments:
            level (int): level of classification desired (between 1 and 4)

        Returns:
            list: list of all EC numbers associated, at the desired level, if
            none are found an empty list is returned
        """
        ec = self.attr.get('EC', None)
        if ec is None:
            return []

        ec = ec.split(',')

        return ['.'.join(x.split('.')[:level]) for x in ec]

    def get_mapping(self, db):
        """
        .. versionadded:: 0.1.13

        Returns the mappings, to a particular db, associated with the
        annotation.

        Arguments:
            db (str): database to which the mappings come from

        Returns:
            list: list of all mappings associated, to the specified db, if
            none are found an empty list is returned
        """
        mappings = self.attr.get('map_{0}'.format(db.upper()))
        if mappings is None:
            return []

        return mappings.split(',')

    def set_mapping(self, db, values):
        """
        .. versionadded:: 0.1.13

        Set mappings to a particular db, associated with the
        annotation.

        Arguments:
            db (str): database to which the mappings come from
            mappings (iterable): iterable of mappings

        """
        self.set_attr(
            'map_{0}'.format(db.upper()),
            ','.join(values)
        )

    @property
    def taxon_id(self):
        """
        taxon_id of the annotation
        """
        value = self.attr.get('taxon_id', None)

        return None if value is None else int(value)

    @taxon_id.setter
    def taxon_id(self, value):
        self.attr['taxon_id'] = int(value)

    @property
    def db(self):
        "db used for the gene_id prediction"
        return self.attr.get('db', None)

    @db.setter
    def db(self, value):
        self.attr['db'] = value

    @property
    def taxon_db(self):
        "db used for the taxon_id prediction"
        return self.attr.get('taxon_db', None)

    @taxon_db.setter
    def taxon_db(self, value):
        self.attr['taxon_db'] = value

    @property
    def dbq(self):
        "db quality of the annotation"
        try:
            return self.get_attr('dbq', int)
        except AttributeNotFound:
            return None

    @dbq.setter
    def dbq(self, value):
        self.attr['dbq'] = value

    @property
    def uid(self):
        """
        .. versionadded:: 0.1.13

        uid of the annotation
        """
        value = self.attr.get('uid', None)
        if value is None:
            #old data where the unique id is marked as ko_idx
            value = self.attr.get('ko_idx', None)

        return value

    @uid.setter
    def uid(self, value):
        self.attr['uid'] = value

    @property
    def bitscore(self):
        "bitscore of the annotation"
        try:
            return float(self.attr['bitscore'])
        except KeyError:
            #legacy for old data
            bitscore = self.attr.get('bit_score', None)
            return None if bitscore is None else float(bitscore)

    @bitscore.setter
    def bitscore(self, value):
        self.attr['bitscore'] = float(value)

    @property
    def gene_id(self):
        "gene_id of the annotation, or *ko* if available"
        try:
            return self.attr['gene_id']
        except KeyError:
            #legacy for old data
            return self.attr.get('ko', None)

    @gene_id.setter
    def gene_id(self, value):
        self.attr['gene_id'] = value

    @property
    def region(self):
        """
        .. versionadded:: 0.1.13

        Return the *region* covered by the annotation, to use in samtools
        """
        return "{0}:{1}:{2}".format(self.seq_id, self.start, self.end)

    def add_exp_syn_count(self, seq, syn_matrix=None):
        """
        .. versionadded:: 0.1.13

        Adds expected synonymous/non-synonymous values for an annotation.

        Arguments:
            seq (str): sequence corresponding to the annotation seq_id
            syn_matrix (None, dict): matrix that determines the return values.
            Defaults to the one defined in the called function
            :func:`mgkit.utils.sequnce.get_seq_expected_syn_count`

        """
        seq = seq[self.start-1:self.end]

        if self.strand == '-':
            seq = seq_utils.reverse_complement(seq)

        syn_count, nonsyn_count = seq_utils.get_seq_expected_syn_count(
            seq,
            syn_matrix=syn_matrix
        )

        self.set_attr('exp_syn', syn_count)
        self.set_attr('exp_nonsyn', nonsyn_count)

    def to_gff(self, sep='='):
        """
        Format the Annotation as a GFF string.

        Arguments:
            sep (str): separator key -> value

        Returns:
            str: annotation formatted as GFF
        """
        var_names = (
            'seq_id', 'source', 'feat_type', 'start', 'end',
            'score', 'strand', 'phase'
        )

        values = '\t'.join(
            str(getattr(self, var_name))
            for var_name in var_names
        )

        attr_column = ';'.join(
            '{0}{1}"{2}"'.format(
                key,
                sep,
                urllib.quote(str(self.attr[key]), ' ()/')
            )
            for key in sorted(self.attr)
        )

        return "{0}\t{1}\n".format(values, attr_column)

    def to_file(self, file_handle):
        """
        Writes the GFF annotation to *file_handle*
        """
        file_handle.write(self.to_gff())

    def to_gtf(self):
        pass

    @property
    def sample_coverage(self):
        """
        .. versionadded:: 0.1.13

        Returns a dictionary with the coverage for each sample, the returned
        dictionary has the sample id (stripped of the *_cov*) suffix and as
        values the coverage (converted via :func:`int`).

        :return dict: dictionary with the samples' coverage
        """
        attributes = self.attr

        return dict(
            (attribute.replace('_cov', ''), int(value))
            for attribute, value in attributes.iteritems()
            if attribute.endswith('_cov')
        )

    def get_number_of_samples(self, min_cov=MIN_COV):
        """
        .. versionadded:: 0.1.13

        Returns the number of sample that have at least a minimum coverage of
        `min_cov`.

        :param int min_cov: minimum coverage
        :return int: number of samples passing the filter
        :raise AttributeNotFound: if no sample coverage attribute is found
        """
        coverage = self.sample_coverage

        if not coverage:
            raise AttributeNotFound(
                'No coverage attribute found (ending in "_cov")'
            )

        return sum(
            1 for sample, coverage in coverage.iteritems()
            if coverage >= min_cov
        )

    def get_attr(self, attr, conv=str):
        """
        .. versionadded:: 0.1.13

        Generic method to get an attribute and convert it to a specific datatype
        """
        try:
            value = self.attr[attr]
        except KeyError:
            raise AttributeNotFound('No {0} attribute found'.format(attr))

        return conv(value)

    def set_attr(self, attr, value):
        """
        .. versionadded:: 0.1.13

        Generic method to set an attribute
        """
        self.attr[attr] = value

    @property
    def coverage(self):
        """
        .. versionadded:: 0.1.13

        Return the total coverage for the annotation

        :return int: coverage
        :raise AttributeNotFound: if no coverage attribute is found
        """
        return self.get_attr('cov', int)

    @property
    def exp_syn(self):
        """
        .. versionadded:: 0.1.13

        Returns the expected number of synonymous changes
        """
        return self.get_attr('exp_syn', int)

    @property
    def exp_nonsyn(self):
        """
        .. versionadded:: 0.1.13

        Returns the expected number of non-synonymous changes
        """
        return self.get_attr('exp_nonsyn', int)

    def get_nuc_seq(self, seq, reverse=False):
        """
        .. versionadded:: 0.1.13

        Returns the nucleotidic sequence that the annotation covers. if the
        annotation's strand is *-*, and *reverse* is True, the reverse
        complement is returned.
        """
        ann_seq = seq[self.start-1:self.end]

        if (self.strand == '-') and reverse:
            ann_seq = seq_utils.reverse_complement(ann_seq)

        return ann_seq


def from_glimmer3(header, line, feat_type='CDS'):
    """
    .. versionadded:: 0.1.12

    Parses the line of a GLIMMER3 ouput and returns an instance of a GFF
    annotation.

    Arguments:
        header (str): the seq_id to which the ORF belongs
        line (str): the prediction line for the orf
        feat_type (str): the feature type to use

    Returns:
        Annotation: instance of annotation

    Example:
        Assuming a GLIMMER3 output like this::

            >sequence0001
            orf00001       66      611  +3     6.08

        The code used is:

        >>> header = 'sequence0001'
        >>> line = 'orf00001       66      611  +3     6.08'
        >>> from_glimmer3(header, line)

    """
    orf_id, start, end, frame, score = line.split()

    start = int(start)
    end = int(end)

    if start > end:
        start, end = end, start

    annotation = Annotation(
        seq_id=header,
        source='GLIMMER3',
        feat_type=feat_type,
        start=start,
        end=end,
        score=float(score),
        strand=frame[0],
        phase=int(frame[1]) - 1,
        frame=frame,
        glimmer_score=float(score),
        orf_id=orf_id
    )

    return annotation


class DuplicateKeyError(Exception):
    """
    .. versionadded:: 0.1.12

    Raised if a GFF annotation contains duplicate keys
    """
    pass


def from_gff(line):
    """
    .. versionadded:: 0.1.12

    Parse GFF line and returns an :class:`Annotation` instance

    Arguments:
        line (str): GFF line

    Returns:
        Annotation: instance of :class:`Annotation` for the line

    Raises:
        DuplicateKeyError: if the attribute column has duplicate keys

    """
    line = line.rstrip()
    line = line.split('\t')

    #in case the last column (attributes) is empty
    if len(line) < 9:
        values = line
        #bug in which the phase was not written
        if len(line[-1]) > 1:
            line.insert(-1, 0)
    else:
        values = line[:-1]

    var_names = (
        'seq_id', 'source', 'feat_type', 'start', 'end',
        'score', 'strand', 'phase'
    )
    #the phase sometimes can be set as unknown, using '-'. We prefer using 0
    var_types = (str, str, str, int, int, float, str, lambda x: 0 if x == '' else int(x))

    attr = {}

    for var, value, vtype in zip(var_names, values, var_types):
        attr[var] = vtype(value)

    #in case the last column (attributes) is empty
    if len(line) < 9:
        return Annotation(**attr)

    for pair in line[-1].split(';'):
        try:
            #by default the key,value separator '=' is assumed to be used
            var, value = pair.strip().split('=', 1)
        except ValueError:
            #in case it doesn't work, it is assumed to be a space
            var, value = pair.strip().split(' ', 1)

        if var in attr:
            raise DuplicateKeyError("Duplicate attribute: {0}".format(var))

        attr[var] = urllib.unquote(value.replace('"', ''))

    return Annotation(**attr)


def from_sequence(name, seq, feat_type='CDS', **kwd):
    """
    .. versionadded:: 0.1.12

    Returns an instance of :class:`Annotation` for the full length of a sequence

    Arguments:
        name (str): name of the sequence
        seq (str): sequence, to get the length of the annotation

    Keyword Args:
        feat_type (str): feature type in the GFF
        **kwd: any additional column

    Returns:
        Annotation: instance of :class:`Annotation`

    """
    annotation = Annotation(
        seq_id=name,
        source='SEQUENCE',
        feat_type=feat_type,
        start=1,
        end=len(seq),
        score=0.0,
        strand='+',
        phase=0,
        sequence=name,
        **kwd
    )

    return annotation


def from_aa_blast_frag(hit, parent_ann, aa_seqs):
    frag_id, frame = hit[0].split('-')
    strand = '+' if frame.startswith('f') else '-'
    frame = int(frame[1])
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]
    if strand == '-':
        start, end = seq_utils.reverse_aa_coord(
            start,
            end,
            len(aa_seqs[hit[0]])
        )
    start, end = seq_utils.convert_aa_to_nuc_coord(start, end, frame)

    annotation = Annotation(
        seq_id=parent_ann.seq_id,
        source='BLAST',
        feat_type='CDS',
        start=start+parent_ann.start-1,
        end=end+parent_ann.start-1,
        score=bitscore,
        strand=strand,
        phase=frame,
        db='UNIPROT',
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        ID=frag_id
    )

    return annotation


def from_nuc_blast_frag(hit, parent_ann, db='NCBI-NT'):
    frag_id = hit[0]
    strand = '+'
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]

    annotation = Annotation(
        seq_id=parent_ann.seq_id,
        source='BLAST',
        feat_type='CDS',
        start=start+parent_ann.start-1,
        end=end+parent_ann.start-1,
        score=bitscore,
        strand=strand,
        phase=0,
        db=db,
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        ID=frag_id
    )

    return annotation


def annotate_sequence(name, seq, window=None):

    length = len(seq)

    if window is None:
        window = length

    for index in xrange(1, length, window):
        annotation = Annotation.from_sequence(name, seq)
        annotation.start = index
        annotation.end = index + window - 1
        if annotation.end > length:
            annotation.end = length
        yield annotation


def from_nuc_blast(hit, db, feat_type='CDS', seq_len=None, **kwd):
    """
    .. versionadded:: 0.1.12

    Returns an instance of :class:`Annotation`

    Arguments:
        hit (tuple): a BLAST hit, from :func:`mgkit.io.blast.parse_blast_tab`
        db (str): db used with BLAST

    Keyword Args:
        feat_type (str): feature type in the GFF
        seq_len (int): sequence length, if supplied, the phase for strand '-'
            can be assigned, otherwise is assigned a 0
        **kwd: any additional column

    Returns:
        Annotation: instance of :class:`Annotation`
    """
    seq_id = hit[0]
    strand = '+'
    identity = hit[2]
    bitscore = hit[-1]
    start = hit[3]
    end = hit[4]

    if start > end:
        start, end = end, start
        strand = '-'
        if seq_len is None:
            phase = 0
        else:
            if (seq_len - end + 1) % 2 == 0:
                phase = 1
            elif (seq_len - end + 1) % 3 == 0:
                phase = 2
            else:
                phase = 0

    if strand == '+':
        if start % 2 == 0:
            phase = 2
        elif start % 3 == 0:
            phase = 2
        else:
            phase = 0

    annotation = Annotation(
        seq_id=seq_id,
        source='BLAST',
        feat_type=feat_type,
        start=start,
        end=end,
        score=bitscore,
        strand=strand,
        phase=phase,
        db=db,
        gene_id=hit[1],
        identity=identity,
        bitscore=bitscore,
        **kwd
    )

    return annotation


def parse_gff(file_handle, gff_type=from_gff):
    """
    .. versionchanged:: 0.1.12
        added *gff_type* parameter

    Parse a GFF file and returns generator of :class:`GFFKegg` instances

    Accepts a file handle or a string with the file name

    Arguments:
        file_handle (str, file): file name or file handle to read from
        gff_type (class): class/function used to parse a GFF annotation

    Yields:
        Annotation: an iterator of :class:`Annotation` instances
    """
    if isinstance(file_handle, str):
        file_handle = mgkit.io.open_file(file_handle, 'r')

    LOG.info(
        "Loading GFF from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    for line in file_handle:
        annotation = gff_type(line)
        yield annotation


def diff_gff(files, key_func=None):
    """
    .. versionadded:: 0.1.12

    Returns a simple diff made between a list of gff files. The annotations are
    grouped using *key_func*, so it depends on it to find similar annotations.

    Arguments:
        files (iterable): an iterable of file handles, pointing to GFF files
        key_func (func): function used to group annotations, defaults to this
            key: *(x.seq_id, x.strand, x.start, x.end, x.gene_id, x.bitscore)*

    Returns:
        dict: the returned dictionary keys are determined by key_func and as
        values lists. The lists elements are tuple whose first element is the
        index of the file, relative to *files* and the second element is the
        line number in which the annotation is. Can be used with the
        :mod:`linecache` module.
    """
    if isinstance(files, str) or len(files) == 1:
        return

    if key_func is None:
        key_func = lambda x: (x.seq_id, x.strand, x.start, x.end, x.gene_id, x.bitscore)

    gff_diff = {}

    for index, file_handle in enumerate(files):
        for lineno, annotation in enumerate(parse_gff(file_handle)):
            key = key_func(annotation)
            try:
                gff_diff[key].append((index, lineno))
            except KeyError:
                gff_diff[key] = [(index, lineno)]

    return gff_diff


def annotation_elongation(ann1, annotations):
    """
    .. versionadded:: 0.1.12

    Given an :class:`Annotation` instance and a list of the instances of the
    same class, returns the longest overlapping range that can be found and the
    annotations that are included in it.

    .. warning::

        annotations are not checked for seq_id and strand

    Arguments:
        ann1 (Annotation): annotation to elongate
        annotations (iterable): iterable of :class:`Annotation` instances

    Returns:
        tuple: the first element is the longest range found, while the the
        second element is a set with the annotations used

    """
    used = set([ann1])

    union = (ann1.start, ann1.end)

    for ann2 in annotations:
        new_union = union_range(union[0], union[1], ann2.start, ann2.end)
        if new_union is not None:
            used.add(ann2)
            union = new_union

    return union, used


def elongate_annotations(annotations):
    """
    .. versionadded:: 0.1.12

    Given an iterable of :class:`Annotation` instances, tries to find the all
    possible longest ranges and returns them.

    .. warning::

        annotations are not checked for seq_id and strand

    Arguments:
        annotations (iterable): iterable of :class:`Annotation` instances

    Returns:
        set: set with the all ranges found
    """

    annotations = sorted(annotations, key=lambda x: x.start)

    ranges = set()

    while len(annotations) > 0:
        ann1 = annotations.pop(0)
        union, used = annotation_elongation(ann1, annotations)
        if union is None:
            ranges.add((ann1.start, ann1.end))
        else:
            annotations = sorted(set(annotations) - used, key=lambda x: x.start)
            ranges.add(union)

    return ranges


def annotation_coverage(annotations, seqs, strand=True):
    """
    .. versionadded:: 0.1.12

    Given a list of annotations and a dictionary where the keys are the sequence
    names referred in the annotations and the values are the sequences
    themselves, returns a number which indicated how much the sequence length is
    "covered" in annotations. If *strand* is True the coverage is strand
    specific.

    Arguments:
        annotations (iterable): iterable of :class:`Annotation` instances
        seqs (dict): dictionary in which the keys are the sequence names and the
            the values are the sequneces
        strand (bool): if True, the values are strand specific (the annotations)
            are grouped by (seq_id, strand) instead of seq_id

    Yields:
        tuple: the first element is the key, (seq_id, strand) if *strand* is
        True or seq_id if *strand* is False, and the coverage is the second
        value.
    """

    if strand:
        key_func = lambda x: (x.seq_id, x.strand)
    else:
        key_func = lambda x: x.seq_id

    annotations = group_annotations(
        annotations,
        key_func=key_func
    )

    for key, key_ann in annotations.iteritems():
        if isinstance(key, str):
            seq_len = len(seqs[key])
        else:
            seq_len = len(seqs[key[0]])

        covered = ranges_length(elongate_annotations(key_ann))

        yield key, covered / seq_len * 100


def group_annotations(annotations, key_func=lambda x: (x.seq_id, x.strand)):
    """
    .. versionadded:: 0.1.12

    Group :class:`Annotation` instances in a dictionary by using a key function
    that returns the key to be used in the dictionary.

    Arguments:
        annotations (iterable): iterable with :class:`Annotation` instances
        key_func (func): function used to extract the key used in the
            dictionary, defaults to a function that returns
            (ann.seq_id, ann.strand)

    Returns:
        dict: dictionary whose keys are returned by *key_func* and the values
        are lists of annotations

    Example:
        >>> ann = [Annotation(seq_id='seq1', strand='+', start=10, end=15),
        ... Annotation(seq_id='seq1', strand='+', start=1, end=5),
        ... Annotation(seq_id='seq1', strand='-', start=30, end=100)]
        >>> group_annotations(ann)
        {('seq1', '+'): [seq1(+):10-15, seq1(+):1-5], ('seq1', '-'): [seq1(-):30-100]}
    """
    grouped = {}

    for annotation in annotations:
        key = key_func(annotation)
        try:
            grouped[key].append(annotation)
        except KeyError:
            grouped[key] = [annotation]

    return grouped


def group_annotations_sorted(annotations, key_func=lambda x: (x.seq_id, x.strand)):
    """
    .. versionadded:: 0.1.13

    Group :class:`Annotation` instances by using a key function that returns a
    key. Assumes that the annotations are already sorted to return an iterator
    and save memory. One way to sort them is using: `sort -s -k 1,1 -k 7,7` on
    the file.

    Arguments:
        annotations (iterable): iterable with :class:`Annotation` instances
        key_func (func): function used to extract the key used in the
            dictionary, defaults to a function that returns
            (ann.seq_id, ann.strand)

    Yields:
        list: a list of the grouped annotations by *key_func* values

    """
    curr_key = ''
    curr_ann = []

    for annotation in annotations:
        new_key = key_func(annotation)
        if curr_key == new_key:
            curr_ann.append(annotation)
        else:
            if curr_key == '':
                curr_ann.append(annotation)
                curr_key = new_key
            else:
                yield curr_ann
                curr_key = new_key
                curr_ann = [annotation]
    else:
        yield curr_ann


def correct_old_annotations(annotations, taxonomy):
    """
    .. versionadded:: 0.1.13

    Corrects old annotations containing a mix of taxonomic annotations (or none
    at all), plus some misspelled taxa.

    * BLAST assigned ID from `blast_taxon_idx` (a number)
    * Profile assigned ID from `taxon_id` which can be in the forms:

        * `id`
        * `name.id`

    * A taxon *name* in which case uses the provided  taxonomy to find its
      ID and returns the first one matching or `None` if no taxonomy is
      passed.

    The taxon_id attribute is set to the correct one, prefferring the
    blast_taxon_idx, then taxon_id, which can be attached to the taxon name and
    as last resort tries to reverse lookup the taxon name.

    :param taxonomy: taxonomy used to resolve the taxon name
    """

    LOG.debug('Correcting old annotations')

    for annotation in annotations:
        #a taxon id from blast
        if ('blast_taxon_idx' in annotation.attr):
            taxon_id = annotation.attr['blast_taxon_idx']
        #a taxon id from profile
        elif 'taxon_id' in annotation.attr:
            taxon_id = annotation.taxon_id
        #a taxon name is provided
        else:
            #if a taxon_name contains the id
            if len(annotation.attr['taxon'].split('.')) == 2:
                taxon_id = annotation.attr['taxon'].split('.')[1]
            #if a taxon_name DOESN'T contains the id try to reverse
            #it using the taxonomy (if provided), using the first matching ID
            else:
                taxon_name = annotation.attr['taxon'].replace('#', ' ')

                if taxon_name in taxon.MISPELLED_TAXA:
                    taxon_name = taxon.MISPELLED_TAXA[taxon_name]

                taxon_id = taxonomy.find_by_name(taxon_name)[0]

        if taxon_id is not None:
            annotation.taxon_id = int(taxon_id)


def extract_nuc_seqs(annotations, seqs, name_func=lambda x: x.uid,
                     reverse=False):
    """
    .. versionadded:: 0.1.13

    Extract the nucleotidic sequences from a list of annotations. Internally
    uses the method :meth:`Annotation.get_nuc_seq`.

    Arguments:
        annotations (iterable): iterable of :class:`Annotation` instances
        seqs (dict): dictionary with the sequences referenced in the annotations
        name_func (func): function used to extract the sequence name to be used,
            defaults to the uid of the annotation
        reverse (bool): if True the annotations on the *-* strand are reverse
            complemented

    Yields:
        tuple: tuple whose first element is the sequence name and the second is
        the sequence to which the annotation refers.
    """
    for annotation in annotations:
        name = name_func(annotation)
        seq = annotation.get_nuc_seq(seqs[annotation.seq_id], reverse=reverse)

        yield name, seq


def group_annotations_by_ancestor(annotations, ancestors, taxonomy):
    """
    .. versionadded:: 0.1.13

    Group annotations by the ancestors provided.

    Arguments:
        annotations (iterable): annotations to group
        ancestors (iterable): list of ancestors accepted
        taxonomy: taxonomy class

    Returns:
        dict: grouped annotations
    """
    ann_dict = dict((ancestor, []) for ancestor in ancestors)

    unknown = []

    for annotation in annotations:
        anc_found = False
        for ancestor, anc_ids in ancestors.iteritems():
            if taxonomy.is_ancestor(annotation.taxon_id, anc_ids):
                ann_dict[ancestor].append(annotation)
                anc_found = True
                break
        if not anc_found:
            unknown.append(annotation)

    return ann_dict, unknown


def split_gff_file(file_handle, name_mask, num_files=2):
    """
    .. versionadded:: 0.1.14

    Splits a GFF, or a list of them, into a number of files. It is assured that
    annotations for the same sequence are kept in the same file, which is
    useful for cases like filtering, even when the annotations are from
    different GFF files.

    Internally, a structure is kept to check if a sequence ID is already been
    stored to a file, in which case the annotation is written to that file,
    otherwise a random file handles (among the open ones) is chosen.

    Arguments:
        file_handle (str, list): a single or list of file handles (or file
           names), from which the GFF annotations are read
        name_mask (str): a string used as template for the output file names
            on which the function applies :func:`string.format`
        num_files (int): the number of files to split the records

    Example:
        >>> import glob
        >>> files = glob.glob('*.gff')
        >>> name_mask = 'split-file-{0}.gff'
        >>> split_gff_file(files, name_mask, 5)
    """
    if isinstance(file_handle, str):
        file_handle = [file_handle]

    file_handle = itertools.chain(
        *(mgkit.io.open_file(x, 'r') for x in file_handle)
    )

    out_handles = [
        open(name_mask.format(filen), 'w')
        for filen in xrange(num_files)
    ]

    seq_ids = {}

    for line in file_handle:
        seq_id = line.split('\t')[0]
        try:
            out_handle = out_handles[seq_ids[seq_id]]
        except KeyError:
            new_index = random.randint(0, num_files - 1)
            seq_ids[seq_id] = new_index
            out_handle = out_handles[new_index]

        out_handle.write(line)
