"""
Module containing functions related to sequence data

.. note::

    For those functions without a docstring, look at the same with a
    underscore ('_') prepended.

"""
from __future__ import division  # add check to use only on python 2.x
from builtins import range, zip
from future.utils import viewitems, viewvalues
import sys
import itertools
import logging
import random
import collections
import numpy
import pandas
from scipy import stats
from scipy import interpolate
import statsmodels.api as sm

from ..utils.common import between
from .trans_tables import UNIVERSAL
from ..io import fasta
from ._sequence import get_kmers, sliding_window, sequence_signature, \
    signatures_matrix

LOG = logging.getLogger(__name__)

TRANS_TABLE = UNIVERSAL
"""Translation table - Universal genetic code"""

REV_COMP = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}
"Dictionary containing the complement of each nucleotide sequence"


def make_reverse_table(tbl=None):
    """
    Makes table to reverse complement a sequence by :func:`reverse_complement`.
    The table used is the complement for each nucleotide, defaulting to
    :data:`REV_COMP`
    """
    if tbl is None:
        tbl = REV_COMP
    trans_table = [chr(idx) for idx in range(256)]
    for nuc, rev in viewitems(tbl):
        trans_table[ord(nuc)] = rev
    return ''.join(trans_table)


# The maketrans function is not available on Python3
if sys.version_info[0] == 2:
    from string import maketrans
    REV_COMP_ASCII = maketrans('ATCG', 'TAGC')
    REV_COMP_UNICODE = {
        ord(unicode(key)): unicode(value)
        for key, value in viewitems(REV_COMP)
    }
else:
    REV_COMP_ASCII = u''.maketrans(REV_COMP)


def reverse_complement_old(seq, tbl=None):
    """
    Returns the reverse complement of a nucleotide sequence

    :param str seq: nucleotide sequence with uppercase characters
    :param dict tbl: dictionary of complement bases, like :data:`REV_COMP`

    :return str: returns the reverse complement of a nucleotide sequence
    """
    if tbl is None:
        tbl = REV_COMP
    rev_complement = [tbl[x] for x in seq[::-1]]
    return ''.join(rev_complement)


def reverse_complement(seq, tbl=REV_COMP_ASCII):
    """
    Returns the reverse complement of a nucleotide sequence

    :param str seq: nucleotide sequence with uppercase characters
    :param dict tbl: translation table returned by :func:`make_reverse_table`

    :return str: returns the reverse complement of a nucleotide sequence
    """
    # Python 2.7 hack with unicode, make sure futurize doesn't change to str
    # since on Python 3 the first condition won't be true, the second won't be
    # evaluated
    if (sys.version_info[0] == 2) and isinstance(seq, unicode):
        tbl = REV_COMP_UNICODE
    return seq[::-1].translate(tbl)


def translate_sequence(sequence, start=0, tbl=None, reverse=False):
    """
    Translate a nucleotide sequence in an amino acid one.

    :param str sequence: sequence to translate, it's expected to be all caps
    :param int start: 0-based index for the translation to start
    :param dict tbl: dictionary with the translation for each codon
    :param bool reverse: if True, :func:`reverse_complement` will be called and
        the returned sequence translated

    :return str: the translated sequence
    """
    if reverse:
        sequence = reverse_complement(sequence)
    if tbl is None:
        tbl = UNIVERSAL

    trs = []

    for idx in range(start, len(sequence), 3):
        codon = sequence[idx:idx+3]
        if len(codon) < 3:
            break
        try:
            trs.append(tbl[codon])
        except KeyError:
            trs.append('X')
    return ''.join(trs)


def put_gaps_in_nuc_seq(nuc_seq, aa_seq, trim=True):
    """
    Match the gaps in an amino-acid aligned sequence to its original nucleotide
    sequence. If the nucleotide sequence is not a multiple of 3, the trim
    option by default trim those bases from the output.

    :param str nuc_seq: original nucleotide sequence
    :param str aa_seq: aligned amino-acid sequence
    :param bool trim: if True trim last nucleotide(s)

    :return str: gapped nucleotide sequence
    """
    # transform a string into a list of codons, then inserts the gaps
    nuc_alg = [nuc_seq[idx:idx+3] for idx in range(0, len(nuc_seq), 3)]

    for idx, aa_pos in enumerate(aa_seq):
        if aa_pos == '-':
            nuc_alg.insert(idx, '-' * 3)

    nuc_alg = ''.join(nuc_alg)

    if trim and (len(nuc_alg) != len(aa_seq) * 3):
        trim_idx = len(nuc_seq) % 3
        if trim_idx > 0:
            nuc_alg = nuc_alg[:-trim_idx]

    return nuc_alg


def get_syn_matrix(trans_table=None, nuc_list=None):
    """
    Returns a dictionary containing the expected count of synonymous and
    non-synonymous changes that a codon can have if one base is allowed to
    change at a time.

    There are 9 possible changes per codon.

    :param dict trans_table: a tranlation table, defaults to
        :data:`seq_utils.TRANS_TABLE`
    :param iterable nuc_list: a list of nucleotides in which a base can change,
        default to the keys of :data:`seq_utils.REV_COMP`

    :return dict: returns a dictionary in which for each codon a dictionary
        {'syn': 0, 'nonsyn': 0} holds the number of expected changes
    """

    if trans_table is None:
        trans_table = TRANS_TABLE
    if nuc_list is None:
        nuc_list = REV_COMP

    syn_matrix = dict(
        (codon, {'syn': 0, 'nonsyn': 0})
        for codon in trans_table
    )

    for codon1 in trans_table:
        for frame in range(0, 3):
            for nuc in nuc_list:
                codon2 = list(codon1)
                codon2[frame] = nuc
                codon2 = ''.join(codon2)

                # skips if its the same
                if codon1 == codon2:
                    continue

                if trans_table[codon1] == trans_table[codon2]:
                    syn_matrix[codon1]['syn'] += 1
                else:
                    syn_matrix[codon1]['nonsyn'] += 1

    return syn_matrix


def get_syn_matrix_all(trans_table=None):
    """
    Same as :func:`get_syn_matrix` but a codon can change in any of the ones
    included in trans_table.

    There are 63 possible changes per codon.
    """

    if trans_table is None:
        trans_table = TRANS_TABLE

    syn_matrix = dict(
        (codon, {'syn': 0, 'nonsyn': 0})
        for codon in trans_table
    )

    for codon1, codon2 in itertools.permutations(trans_table, 2):
        if trans_table[codon1] == trans_table[codon2]:
            syn_matrix[codon1]['syn'] += 1
        else:
            syn_matrix[codon1]['nonsyn'] += 1

    return syn_matrix


def get_seq_expected_syn_count(seq, start=0, syn_matrix=None):
    """
    Calculate the expected number of synonymous and non-synonymous changes in a
    nucleotide sequence. Assumes that the sequence is already in the correct
    frame and its length is a multiple of 3.

    :param iterable seq: nucleotide sequence (uppercase chars)
    :param int start: frame of the sequence
    :param dict syn_matrix: dictionary that contains the expected number of
        changes for a codon, as returned by :func:`get_syn_matrix`

    :return tuple: tuple with counts of expected counts (syn, nonsyn)
    """
    if syn_matrix is None:
        syn_matrix = get_syn_matrix()

    seq_syn = 0
    seq_nonsyn = 0

    for idx in range(start, len(seq), 3):
        codon = seq[idx:idx+3]
        seq_syn += syn_matrix[codon]['syn']
        seq_nonsyn += syn_matrix[codon]['nonsyn']

    return seq_syn, seq_nonsyn


def get_variant_sequence(seq, *snps):
    """
    .. versionadded:: 0.1.16

    Return a sequence changed in the positions requested.

    Arguments:
        seq (str): a sequence
        *snps (tuple): each argument passed is a tuple with the first element
            as a position in the sequence (1-based index) and the second
            element is the character to substitute in the sequence

    Returns:
        str: string with the changed characters

    Example:
        >>> get_variant_sequence('ACTGATATATGCGCGCATCT', (1, 'C'))
        'CCTGNTGTATGCGCGCATCT'

    .. note::

        It is used for nucleotide sequences, but it is valid to use any string
    """
    seq = list(seq)

    for pos, change in snps:
        seq[pos - 1] = change

    return ''.join(seq)


def convert_aa_to_nuc_coord(start, end, frame=0):
    """
    Converts aa coordinates to nucleotidic ones. The coordinates must be from
    '+' strand. For the '-' strand, use :func:`reverse_aa_coord` first.

    Arguments:
        start (int): start of the annotation (lowest number)
        end (int): end of the annotation (highest number)
        frame (int): frame of the AA translation (0, 1 or 2)

    Returns:
        tuple: the first element is the converted *start* and the second
        element is the converted *end*

    .. note::

        the coordinates are assumed to be 1-based indices

    """
    start = (start - 1) * 3 + 1  # gets the first base of the codon
    end = (end - 1) * 3 + 3  # gets the third base of the codon

    return (start + frame, end + frame)


def reverse_aa_coord(start, end, seq_len):
    """
    Used to reverse amino-acid coordinates when parsing an AA annotation on
    the - strand. Used when the BLAST or HMMER annotations use AA sequences.

    Arguments:
        start (int): start of the annotation
        end (int): end of the annotation
        seq_len (int): aa sequence length

    Returns:
        tuple: reversed (from strand - to strand +) coordinates. The first
        element is the converted *start* and the second element is the
        converted *end*

    .. note::

        * start and end are 1-based indices

    """

    return (seq_len - end + 1, seq_len - start + 1)


def calc_n50(seq_lengths):
    """
    Calculate the N50 statistics for a :class:`numpy.array` of sequence
    lengths.

    The algorithm finds in the supplied array the element (contig length) for
    which the sum all contig lengths equal or greater than it is equal to half
    of all assembled base pairs.

    :param array seq_lengths: an instance of a numpy array containing the
        sequence lengths

    :return int: the N50 statistics value
    """
    LOG.info("Calculating N50 for assembly")
    seq_lengths.sort()

    bp_half = seq_lengths.sum() / 2.0
    bp_sum = 0

    for seq_length in seq_lengths[::-1]:
        bp_sum += seq_length
        if bp_sum > bp_half:
            return seq_length


class Alignment(object):
    "Simple alignment class"
    _names = None
    _seqs = None

    def __init__(self, seqs=None):
        self._names = []
        self._seqs = []
        if seqs:
            self.add_seqs(seqs)

    def add_seq(self, name, seq):
        """
        Add a sequence to the alignment

        :param str name: name of the sequence
        :param str seq: sequence
        """
        self._names.append(name)
        self._seqs.append(seq)

    def add_seqs(self, seqs):
        """
        Add sequences to the alignment

        :param iterable seqs: iterable that returns (name, seq)
        """
        for name, seq in seqs:
            self.add_seq(name, seq)

    def get_position(self, pos):
        """
        Get all characters at a position

        :param int pos: position to return (0-based)

        :return str: all characters occuring at the position
        """
        return ''.join(seq[pos] for seq in self._seqs)

    def get_seq_len(self):
        "Get the length of the alignment"
        return len(self._seqs[0])

    def get_consensus(self, nucl=True):
        """
        .. versionchanged:: 0.1.16
            added *nucl* parameter

        The consensus sequence is constructed by checking the nucleotide that
        has the maximum number of counts for each position in the alignment.

        Arguments:
            nucl (bool): specify if the alignment is nucleotidic

        Returns:
            str: consensus sequence
        """

        cons_seq = []
        nucs = REV_COMP.keys()

        for pos in range(0, self.get_seq_len()):
            site = self.get_position(pos)
            if nucl:
                nuc_pos = max(
                    ((site.count(nuc), nuc) for nuc in nucs),
                    key=lambda x: x[0]
                )[1]
            else:
                nuc_pos = max(collections.Counter(x for x in site if x != '-'))

            cons_seq.append(nuc_pos)
        return ''.join(cons_seq)

    def __iter__(self):

        for name, seq in zip(self._names, self._seqs):
            yield name, seq

    def get_snps(self, ref_seq=None, full_size=False):
        """
        A SNP is called for the nucleotide that has the most counts among the
        ones that differ in the each site of the alignment. If two nucleotides
        have the same maximum count, one is randomly chosen.

        :param str ref_seq: a reference sequence can be provided, if None, a
            consensus sequence is produced for the alignment
        :param bool full_size: if True a tuple is returned for each position in
            the alignment. If there is no SNP at a position the value for the
            SNP is None

        :return list: a list of tuples (position, SNP)
        """
        if ref_seq is None:
            ref_seq = self.get_consensus()

        nucs = REV_COMP.keys()
        snps = []

        for pos in range(0, self.get_seq_len()):
            ref_nuc = ref_seq[pos]

            alts = self.get_position(pos).replace(ref_nuc, '').replace('-', '')
            if not alts:
                if full_size:
                    snps.append((pos, None))
                continue
            alts = sorted(
                ((alts.count(nuc), nuc) for nuc in nucs),
                key=lambda x: x[0]
            )
            counts = [alt[0] for alt in alts]
            count_idx = counts.index(counts[-1])
            alt = random.choice(alts[count_idx:])
            snps.append((pos, alt[1]))

        return snps


def get_seq_number_of_syn(ref_seq, snps, start=0, trans_table=None):
    """
    Given a reference sequence and a list of SNPs, calculates the number of
    synonymous and non-synonymous SNP.

    :param str ref_seq: reference sequence
    :param iterable snps: list of tuples (position, SNP) - zero based index
    :param int start: the frame used for the reference {0, 1, 2}
    :param dict trans_table: translation table used - codon->AA

    :return tuple: synonymous and non-synonymous counts
    """
    if trans_table is None:
        trans_table = TRANS_TABLE

    syn = 0
    nonsyn = 0

    for idx in range(start, len(ref_seq), 3):
        ref_codon = ref_seq[idx:idx+3]
        end = idx + 2
        for pos, change in snps:
            if not between(pos, idx, end):
                continue
            change_pos = pos-idx
            codon = ref_codon[:change_pos] + change + ref_codon[change_pos+1:]

            if len(codon) < 3:
                continue

            if trans_table[codon] == trans_table[ref_codon]:
                syn += 1
            else:
                nonsyn += 1

    return syn, nonsyn


def check_snp_in_seq(ref_seq, pos, change, start=0, trans_table=None):
    """
    Check a SNP in a reference sequence if it is a synonymous or non-synonymous
    change.

    :param str ref_seq: reference sequence
    :param int pos: SNP position - it is expected to be a 1 based index
    :param str change: nucleotide change occuring at *pos*
    :param int start: the starting position for the coding region - 0 based
        index
    :param dict trans_table: translation table used - codon->AA

    :return bool: True if it is a synonymous change, False if non-synonymous
    """
    # pos is expected to be a 1 based index
    # returns true if snp is synonymous
    # start is expected to be 0 based
    if trans_table is None:
        trans_table = TRANS_TABLE

    codon_idx = ((pos - start) // 3) - (1 if (pos - start) % 3 == 0 else 0)

    codon_pos = codon_idx * 3 + start

    ref_codon = ref_seq[codon_pos:codon_pos + 3]
    snp_codon = list(ref_codon)
    snp_codon[(pos % 3) - 1] = change
    snp_codon = ''.join(snp_codon)

    return trans_table[snp_codon] == trans_table[ref_codon]


def sequence_gc_ratio(sequence):
    """
    .. versionadded:: 0.1.13

    Calculate GC ratio information for a sequence. The formula is:

    .. math::
        :label: gc_ratio

        \\frac {(A + T)}{(G + C)}

    Arguments:
        sequence (str): sequence

    Returns:
        float: GC ratio, or `numpy.nan` if G = C = 0
    """

    at_sum = (sequence.count('A') + sequence.count('T'))
    gc_sum = (sequence.count('G') + sequence.count('C'))

    try:
        return at_sum / gc_sum
    except ZeroDivisionError:
        return numpy.nan


def sequence_gc_content(sequence):
    """
    .. versionchanged:: 0.3.3
        in case of `ZeroDivisionError` returns .5

    .. versionadded:: 0.1.13

    Calculate GC content information for an annotation. The formula is:

    .. math::
        :label: gc_content

        \\frac {(G + C)}{(G + C + A + T)}


    Arguments:
        sequence (str): sequence

    Returns:
        float: GC content
    """

    at_sum = (sequence.count('A') + sequence.count('T'))
    gc_sum = (sequence.count('G') + sequence.count('C'))

    try:
        return gc_sum / (gc_sum + at_sum)
    except ZeroDivisionError:
        return .5


def sequence_composition(sequence, chars=tuple(REV_COMP)):
    """
    .. versionadded:: 0.1.13

    Returns the number of occurences of each unique character in the sequence

    Arguments:
        sequence (str): sequence
        chars (iterable, None): iterable of the chars to test, default to
            (A, C, T, G). if None checks all unique characters in the sequencce

    Yields:
        tuple: the first element is the nucleotide and the second is the number
        of occurences in *sequence*
    """
    if chars is None:
        chars = set(sequence)

    for nuc in chars:
        yield nuc, sequence.count(nuc)


def get_contigs_info(file_name, pp=False):
    """
    .. versionchanged:: 0.2.4
       file_name can be a *dict* name->seq or a list of sequences

    .. versionadded:: 0.2.1

    Given a file name for a fasta file with sequences, a dictionary of
    name->seq, or a list of sequences, returns the following information in a
    tuple, or a string if *pp* is True:

    - number of sequences
    - total base pairs
    - max length
    - min length
    - average length
    - N50 statistic

    Arguments:
        file_name (str): fasta file to open
        pp (bool): if True, a formatted string is returned

    Returns:
        str, tuple: the returned value depends on the value of *pp*, if True a
        formatted string is returned, otherwise the tuple with all values is.
    """

    if isinstance(file_name, dict):
        seqs = list(viewvalues(file_name))
        file_name = 'dictionary'
    elif isinstance(file_name, list):
        seqs = list(file_name)
        file_name = 'list'
    else:
        seqs = list(seq for name, seq in fasta.load_fasta(file_name))

    lengths = numpy.array([len(x) for x in seqs])

    info = (
        len(seqs),
        lengths.sum(),
        lengths.max(),
        lengths.min(),
        lengths.mean(),
        calc_n50(lengths),
    )

    if pp:
        info_str = '{}: {} contigs, {} bp, max {}, min {}, avg {:.2f}, n50 {}'
        info = info_str.format(
            file_name, *info
        )

    return info


def _sliding_window(seq, size, step=None):
    """
    .. versionadded:: 0.2.6

    Returns a generator, with every iteration yielding a subsequence of size
    *size*, with a step of *step*.

    Arguments:
        seq (str): sequnece
        size (int): size of the sliding window
        step (int, None): the step to use in the sliding window. If *None*,
            half of the sequence length is used

    Yields:
        str: a subsequence of size *size* and step *step*
    """
    for index in range(0, len(seq) - size + 1, size // 2 if step is None else step):
        yield seq[index:index+size]


def _get_kmers(seq, k):
    """
    .. versionadded:: 0.2.6

    Returns a generator, with every iteration yielding a kmer of size *k*

    Arguments:
        seq (str): sequence
        k (int): kmer size

    Yields:
        str: a portion of *seq*, of size *k* with a step of *1*
    """
    for index in range(0, len(seq) - k + 1):
        yield seq[index:index+k]


def _sequence_signature(seq, w_size, k_size=4, step=None):
    """
    .. versionadded:: 0.2.6

    Returns the signature of a sequence, based on a kmer length, over a sliding
    window. Each sliding window signature is placed in order into a list, with
    each element being a :class:`collections.Counter` instance whose keys are
    the kmer found in that window.

    Arguments:
        seq (str): sequence for which to get the signature
        w_size (int): size of the sliding window size
        k_size (int): size of the kmer to use :func:`get_kmers`
        step (int): step to use in :func:`sliding_window`

    Returns:
        list: a list of :class:`collections.Counter` instances, for each
        window used
    """
    kmer_counts = []
    for subseq in _sliding_window(seq, w_size, step):
        kmer_counts.append(
            collections.Counter(kmer for kmer in _get_kmers(subseq, k_size) if 'N' not in kmer)
        )
    return kmer_counts


def _signatures_matrix(seqs, w_size, k_size=4, step=None):
    """
    .. versionadded:: 0.2.6

    Return a matrix (pandas.DataFrame) where the columns are the kmer found in
    all sequences *seqs* and the rows are the a MultiIndex with the first level
    being the sequnce name and the second the index of the sliding window for
    which a signature was computed.

    Arguments:
        seqs (iterable): iterable that yields a tuple, with the first element
            being the sequence name and the second the sequence itself
        w_size (int): size of the sliding window size
        k_size (int): size of the kmer to use :func:`get_kmers`
        step (int): step to use in :func:`sliding_window`, defaults to half of
            the window size

    Returns:
        pandas.DataFrame: a DataFrame where the columns are the kmers and the
        rows are the signatures of each contigs/windows.
    """

    def flatten_contigs(data):
        for name, windows in viewitems(data):
            for index, window in enumerate(windows):
                yield (name, index), dict(window)

    step = w_size // 2 if step is None else step

    kmer_counts = {}

    for name, seq in seqs:
        sign = sequence_signature(seq, w_size, k_size, step=step)
        kmer_counts[name] = sign

    return pandas.DataFrame(
        dict(
            flatten_contigs(kmer_counts)
        )
    ).T.fillna(0)


def random_sequences_codon(n=1, length=150, codons=list(UNIVERSAL.keys()),
                           p=None, frame=None):
    """
    .. versionadded:: 0.3.3

    Returns an iterator of nucleotidic sequences, based on a defined genetic
    code (passed as parameter, defaults to the universal one). The sequence if
    first sampled with replacement from the codon list, with a number of codons
    that covers the length chosen plus an additional one to allow a frame shift
    as set by *frame*

    .. note::

        If the probability (for each codon) are supplied, the number of
        sequences required to match those probabilities within a 10% margin of
        error is of at least 10.000 sequences, for 5% at leas 100.000

    Arguments:
        n (int): number of sequences to yield
        length (int): length of the sequences
        codons (iterable): codons used when generating the sequences
        p (tuple): probability of each codon occurence, in the same order as
            *codons*
        frame (int or None): used to define a specific frame shift occuring in
            the sequence (0 to 2) or a random one (if *None*)

    Yields:
        str: string representing a nucleotidic sequence
    """
    sample_size = (length // 3) + (1 if length % 3 != 0 else 0) + 1
    codons = numpy.array(codons)

    if frame is not None:
        if (frame < 0) or (frame > 2):
            frame = 0
        sframe = frame

    for i in range(n):

        if frame is None:
            sframe = numpy.random.randint(0, 4)

        yield ''.join(
            numpy.random.choice(codons, size=sample_size, replace=True, p=p)
        )[sframe:length+sframe]


def random_sequences(n=1, length=150, p=None):
    """
    .. versionadded:: 0.3.3

    Returns an iterator of random squences, where each nucleotide probability
    can be customised in the order (A, C, T, G)

    Arguments:
        n (int): number of sequences to yield
        length (int): length of each sequence
        p (tuple): tuple with the probability of a nucleotide to occur, in the
            order A, C, T, G

    Yields:
        str: string representing a nucleotidic sequence
    """
    nucl = numpy.array(['A', 'C', 'T', 'G'])

    for i in range(n):

        yield ''.join(
            numpy.random.choice(nucl, size=length, replace=True, p=p)
        )


def qualities_model_decrease(length=150, scale=None, loc=35):
    """
    .. versionadded:: 0.3.3

    The model is a decreasing one, from 35 and depends on the length of the
    sequence.

    Arguments:
        length (int): length of the qualities
        scale (float): base level of the qualities
        loc (float): loc parameter of the normal distribution

    Return:
        tuple: first element is sequence qualities, the second element contains
        the distribution used to randomise them
    """
    if scale is None:
        scale = 2. / (length / 150)
    return loc - (numpy.log(numpy.arange(length) + 1) * (length / 150)), stats.norm.freeze(loc=0, scale=scale)


def qualities_model_constant(length=150, scale=1, loc=35):
    """
    .. versionadded:: 0.3.3

    Model with constant quality

    Arguments:
        length (int): length of the qualities
        scale (float): base level of the qualities
        loc (float): loc parameter of the normal distribution

    Return:
        tuple: first element is sequence qualities, the second element contains
        the distribution used to randomise them
    """
    return numpy.ones(length) * loc, stats.norm.freeze(scale=scale)


def extrapolate_model(quals, frac=.5, scale_adj=.5):
    """
    .. versionadded:: 0.3.3

    Extrapolate a quality model from a list of qualities. It uses internally
    a LOWESS as the base, which is used to estimate the noise as a normal
    distribution.

    Arguments:
        quals (list): list of arrays of qualities, sorted by position in the
            corresponding sequence
        frac (float): fraction of the data used for the LOWESS fit (uses
            statsmodels)
        scale_adj (float): value by which the scale of the normal distribution
            will be multiplied. Defaults to halving the scale

    Returns:
        tuple: the first element is the qualities fit with a LOWESS, the second
        element is the distribution
    """
    if not isinstance(quals, list):
        quals = list(quals)

    endog = numpy.hstack(quals)
    exog = numpy.hstack(
        [numpy.arange(len(qual)) + 1 for qual in quals]
    )

    lowess = sm.nonparametric.lowess(endog, exog, frac=frac)
    lowess = interpolate.interp1d(
        lowess[:, 0],
        lowess[:, 1],
        kind='nearest',
        bounds_error=False,
        fill_value='extrapolate'
    )
    lowess = lowess(numpy.arange(exog.max()) + 1)
    dist = numpy.hstack([
        qual - lowess[:len(qual)]
        for qual in quals
    ])
    dist_args = stats.norm.fit(dist)
    dist_args = (dist_args[0], dist_args[1] * scale_adj)
    dist = stats.norm(*dist_args)
    return lowess, dist


def random_qualities(n=1, length=150, model=None):
    """
    .. versionadded:: 0.3.3

    Arguments:
        n (int): number of quality arrays to yield
        length (int): length of the quality array
        model (tuple): a tuple specifying the qualities and error distribution,
            if *None* :func:`qualities_model_decrease` is used

    Yields:
        numpy.array: numpy array of qualities, with the maximum value of 40
    """
    if model is None:
        model = qualities_model_decrease(length=length)

    base, dist = model

    for x in range(n):
        qual = numpy.round(base + dist.rvs(size=length)).astype(int)
        qual[qual > 40] = 40
        yield qual
