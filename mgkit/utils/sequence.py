"""
Module containing functions related to sequence data
"""
from __future__ import division  # add check to use only on python 2.x

import itertools
import logging
import random
import numpy
from ..utils.common import between
from .trans_tables import UNIVERSAL
import collections

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
    for nuc, rev in tbl.iteritems():
        trans_table[ord(nuc)] = rev
    return ''.join(trans_table)


REV_COMP_ASCII = make_reverse_table()


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

    for idx in xrange(start, len(sequence), 3):
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
    sequence. If the nucleotide sequence is not a multiple of 3, the trim option
    by default trim those bases from the output.

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

    syn_matrix = dict((codon, {'syn': 0, 'nonsyn': 0}) for codon in trans_table)

    for codon1 in trans_table:
        for frame in range(0, 3):
            for nuc in nuc_list:
                codon2 = list(codon1)
                codon2[frame] = nuc
                codon2 = ''.join(codon2)

                #skips if its the same
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

    syn_matrix = dict((codon, {'syn': 0, 'nonsyn': 0}) for codon in trans_table)

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
        seq_len (int): aa sequence lenght

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
    Calculate the N50 statistics for a :class:`numpy.array` of sequence lengths.

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
    #pos is expected to be a 1 based index
    #returns true if snp is synonymous
    #start is expected to be 0 based
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

    return gc_sum / (gc_sum + at_sum)


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
