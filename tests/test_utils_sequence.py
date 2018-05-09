from __future__ import division
from builtins import zip
import pytest
import numpy
import itertools
import pandas
from collections import Counter
from mgkit.utils import trans_tables
from mgkit.utils.sequence import reverse_aa_coord, get_variant_sequence, \
    convert_aa_to_nuc_coord, reverse_complement, translate_sequence, \
    put_gaps_in_nuc_seq, get_seq_expected_syn_count, Alignment, \
    get_seq_number_of_syn, sequence_gc_ratio, sequence_gc_content, \
    sequence_composition, _get_kmers, _sliding_window, _sequence_signature, \
    sequence_signature, sliding_window, get_kmers, random_sequences_codon, \
    random_sequences


@pytest.mark.parametrize(
    "coords,result",
    [
        ((6, 8, 8), (1, 3)),
        ((2, 4, 8), (5, 7)),
        ((3, 5, 8), (4, 6)),
        ((1, 5, 8), (4, 8)),
        ((5, 8, 8), (1, 4)),
        ((10, 30, 100), (71, 91))
    ]
)
def test_reverse_aa_coord(coords, result):
    assert reverse_aa_coord(*coords) == result


@pytest.mark.parametrize(
    "seq,snps,result",
    [
        ['ACTGATATATGCGCGCATCT', [(1, 'C')], 'CCTGATATATGCGCGCATCT'],
        ['ACTGATATATGCGCGCATCT', [(1, 'C'), (7, 'G'), (5, 'N')], 'CCTGNTGTATGCGCGCATCT'],
        ['ACTGATATATGCGCGCATCT', [(1, 'C'), (7, 'G'), (5, 'N')], 'CCTGNTGTATGCGCGCATCT']
    ]
)
def test_get_variant_sequence(seq, snps, result):
    assert get_variant_sequence(seq, *snps) == result


@pytest.mark.parametrize(
    "coords,result",
    [
        [(1, 4, 0), (1, 12)],
        [(1, 4, 1), (2, 13)],
        [(3, 12, 0), (7, 36)],
        [(2, 23, 2), (6, 71)],
    ]
)
def test_convert_aa_to_nuc_coord(coords, result):
    assert convert_aa_to_nuc_coord(*coords) == result


def test_reverse_complement():
    seq = 'ACTGATATATGCGCGCATCTC'
    rev = 'GAGATGCGCGCATATATCAGT'

    assert reverse_complement(seq) == rev



@pytest.mark.parametrize(
    "seq,trl,start,tbl,reverse",
    [
        ('TTTAAAACCGGGGTC', 'FKTGV', 0, trans_tables.UNIVERSAL, False),
        ('TTTAAAACCGGGGTC', 'FKTGV', 0, trans_tables.VT_MIT, False),
    ]
)
def test_translate_sequence_eq(seq, trl, start, tbl, reverse):
    assert translate_sequence(seq, start=start, tbl=tbl, reverse=reverse) == trl

@pytest.mark.parametrize(
    "seq,trl,start,tbl,reverse",
    [
        ('TTTAAAACCGGGGTC', 'FKTGV', 1, trans_tables.UNIVERSAL, False),
        ('TTTAAAACCGGGGTC', 'FKTGV', 2, trans_tables.UNIVERSAL, False),
    ]
)
def test_translate_sequence_ne(seq, trl, start, tbl, reverse):
    assert translate_sequence(seq, start=start, tbl=tbl, reverse=reverse) != trl


@pytest.mark.parametrize(
    "nuc,aa_gap,nuc_gap, trim",
    [
        ('TTTAAAACCGGGGTC', 'FK-TG-V', 'TTTAAA---ACCGGG---GTC', True),
        ('TTTAAAACCGGGGTC', '-FKTG-V', '---TTTAAAACCGGG---GTC', True),
        ('TTTAAAACCGGGGTC', 'FKTG-V-', 'TTTAAAACCGGG---GTC---', True),
        ('TTTAAAACCGGGGTC', '-FKTG-V-', '---TTTAAAACCGGG---GTC---', True),
        # trim sequence by default
        ('TTTAAAACCGGGGTCAA', '-FKTG-V-', '---TTTAAAACCGGG---GTC---', True),
        # don't trim sequence by default
        ('TTTAAAACCGGGGTCAA', '-FKTG-V-', '---TTTAAAACCGGG---GTC---AA', False),
    ]
)
def test_put_gaps_in_nuc_seq(nuc, aa_gap, nuc_gap, trim):
    assert put_gaps_in_nuc_seq(nuc, aa_gap, trim) == nuc_gap


def test_get_seq_expected_syn_count():
    seq = 'ATGCATCGACTCTGCACTACG'

    assert get_seq_expected_syn_count(seq) == (15, 48)


@pytest.fixture
def alignment():
    return Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCCCT')
        ]
    )


def test_alignment1(alignment):
    assert alignment.get_position(6) == 'TCC'


def test_alignment2(alignment):
    assert alignment.get_seq_len() == 9


def test_alignment3(alignment):
    assert alignment.get_consensus() == 'ACCCACCCT'


def test_alignment4(alignment):
    ref_seq = alignment.get_consensus()
    snps = [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')]
    assert alignment.get_snps(ref_seq) == snps


def test_alignment5(alignment):
    ref_seq = alignment.get_consensus()
    snps = [
        (0, None),
        (1, None),
        (2, 'T'),
        (3, None),
        (4, 'G'),
        (5, None),
        (6, 'T'),
        (7, None),
        (8, 'A')
    ]
    assert alignment.get_snps(ref_seq, full_size=True) == snps


@pytest.mark.parametrize(
    "ref_seq,snps,start,trans_table,result",
    [
        ('ACCCACCCT', [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')], 0, trans_tables.UNIVERSAL, (2, 2)),
        ('ACCCACCCT', [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')], 1, trans_tables.UNIVERSAL, (1, 2)),
        ('ACCCACCCT', [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')], 2, trans_tables.UNIVERSAL, (1, 2)),
    ]
)
def test_get_seq_number_of_syn(ref_seq, snps, start, trans_table, result):
    assert get_seq_number_of_syn(ref_seq, snps, start=start, trans_table=trans_table) == result


def test_sequence_gc_ratio1():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 4 + 'G' * 11
    assert sequence_gc_ratio(seq) == 14 / 15


def test_sequence_gc_ratio_fail1():
    seq = 'A' * 10 + 'T' * 4
    assert numpy.isnan(sequence_gc_ratio(seq))


def test_sequence_gc_content1():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 4 + 'G' * 11
    sequence_gc_content(seq) == 15 / 29


def test_sequence_composition1():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 5 + 'G' * 11
    assert sorted(sequence_composition(seq)) == [('A', 10), ('C', 4), ('G', 11), ('T', 5)]


def test_sequence_composition2():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 5 + 'G' * 11 + 'N' * 2
    assert sorted(sequence_composition(seq, chars=None)) == [('A', 10), ('C', 4), ('G', 11), ('N', 2), ('T', 5)]


def test_sequence_composition3():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 5 + 'G' * 11 + 'N' * 2
    assert sorted(sequence_composition(seq, chars=('A', 'C'))) == \
        [('A', 10), ('C', 4)]


def test_get_kmers1():
    seq = 'ACTG' * 2
    assert list(_get_kmers(seq, 4)) == ['ACTG', 'CTGA', 'TGAC', 'GACT', 'ACTG']


def test_get_kmers2():
    seq = 'ACTG' * 2
    assert list(_get_kmers(seq, 5)) == ['ACTGA', 'CTGAC', 'TGACT', 'GACTG']


def test_sliding_window1():
    seq = 'ACTG' * 5
    assert list(_sliding_window(seq, 4, 4)) == ['ACTG'] * 5


def test_sliding_window2():
    seq = 'ACTG' * 2
    assert list(_sliding_window(seq, 4, 2)) == ['ACTG', 'TGAC', 'ACTG']


def test_sliding_window3():
    seq = 'ACTG' * 2
    assert list(_sliding_window(seq, 4, 3)) == ['ACTG', 'GACT']


def test_sequence_signature1():
    seq = 'ACTG' * 2
    count = _sequence_signature(seq, 4, 4, 4)
    assert len(count) == 2


def test_sequence_signature2():
    seq = 'ACTG' * 2
    count = _sequence_signature(seq, 4, 4, 4)
    assert count[0]['ACTG'] == 1


def test_sequence_signature3():
    seq = 'ACTG' * 2
    count = _sequence_signature(seq, 8, 4, 4)
    assert count[0]['ACTG'] == 2


def test_sequence_signature4():
    seq = 'ACTG' * 2
    count = _sequence_signature(seq, 5, 4, 5)
    assert len(count) == 1


def test_sequence_signature_cython():
    seq = 'ACTG' * 2
    countP = _sequence_signature(seq, 4, 4, 4)
    countC = sequence_signature(seq, 4, 4, 4)
    assert len(countP) == len(countC)


def test_sliding_window_cython():
    seq = 'ACTG' * 5
    assert list(
         _sliding_window(seq, 4, 4)
    ) == list(
         sliding_window(seq, 4, 4)
    )


def test_get_kmers_cython():
    seq = 'ACTG' * 2
    assert list(_get_kmers(seq, 4)) == list(get_kmers(seq, 4))


def test_random_sequences_codon1():
    seqs = random_sequences_codon(n=20, length=150, frame=0, codons=list(trans_tables.UNIVERSAL))
    it = itertools.chain(
        *(
            _sliding_window(seq, 3, 3)
            for seq in seqs
        )
    )
    p = pandas.Series(Counter(it))
    p = p.mean() / p.sum()

    exp_p = 1. / len(trans_tables.UNIVERSAL)
    assert p == exp_p


def test_random_sequences_codon2():
    seqs = random_sequences_codon(n=20, length=150, frame=0, codons=list(trans_tables.VT_MIT))
    it = itertools.chain(
        *(
            _sliding_window(seq, 3, 3)
            for seq in seqs
        )
    )
    p = pandas.Series(Counter(it))
    p = p.mean() / p.sum()

    exp_p = 1. / len(trans_tables.VT_MIT)
    assert p == exp_p


def test_random_sequences_codon3():
    exp_p = pandas.Series(
        numpy.random.normal(
            size=len(trans_tables.UNIVERSAL), loc=100, scale=50
        ),
        index=list(trans_tables.UNIVERSAL.keys())
    ).apply(abs)
    exp_p = (exp_p / exp_p.sum()).sort_index()

    seqs = random_sequences_codon(n=10**5, length=150, frame=0,
                                  codons=sorted(trans_tables.UNIVERSAL),
                                  p=exp_p.tolist())
    it = itertools.chain(
        *(
            _sliding_window(seq, 3, 3)
            for seq in seqs
        )
    )
    p = pandas.Series(Counter(it), index=exp_p.index)
    p = p / p.sum()

    assert all(p < (exp_p * 1.1)) and all(p > (exp_p * .9))
