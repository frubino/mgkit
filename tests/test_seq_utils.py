from nose.tools import *

from mgkit.utils import sequence
from mgkit.utils import trans_tables


def test_reverse_aa_coord():
    res = []
    exp = [(1, 3), (5, 7), (4, 6), (4, 8), (1, 4)]

    res.append(sequence.reverse_aa_coord(6, 8, 8))
    res.append(sequence.reverse_aa_coord(2, 4, 8))
    res.append(sequence.reverse_aa_coord(3, 5, 8))
    res.append(sequence.reverse_aa_coord(1, 5, 8))
    res.append(sequence.reverse_aa_coord(5, 8, 8))

    eq_(res, exp)


def test_convert_aa_to_nuc_coord():
    res = []
    exp = [(1, 12), (2, 13), (7, 36), (6, 71)]

    res.append(sequence.convert_aa_to_nuc_coord(1, 4, 0))
    res.append(sequence.convert_aa_to_nuc_coord(1, 4, 1))

    res.append(sequence.convert_aa_to_nuc_coord(3, 12, 0))
    res.append(sequence.convert_aa_to_nuc_coord(2, 23, 2))

    eq_(res, exp)


def test_reverse_complement():
    seq = 'ACTGATATATGCGCGCATCT'
    rev = 'AGATGCGCGCATATATCAGT'

    eq_(rev, sequence.reverse_complement(seq))


def test_translate_sequence1():
    seq = 'TTTAAAACCGGGGTC'
    trl = 'FKTGV'

    tseq = sequence.translate_sequence(seq, start=0, tbl=trans_tables.UNIVERSAL)

    eq_(trl, tseq)


def test_translate_sequence2():
    seq = 'TTTAAAACCGGGGTC'
    trl = 'FKTGV'

    tseq = sequence.translate_sequence(seq, start=1, tbl=trans_tables.UNIVERSAL)

    ok_(trl != tseq)


def test_translate_sequence3():
    seq = 'TTTAAAACCGGGGTC'
    trl = 'FKTGV'

    tseq = sequence.translate_sequence(seq, start=2, tbl=trans_tables.UNIVERSAL)

    ok_(trl != tseq)


def test_put_gaps_in_nuc_seq1():
    seq = 'TTTAAAACCGGGGTC'
    trl = 'FK-TG-V'

    eq_(sequence.put_gaps_in_nuc_seq(seq, trl), 'TTTAAA---ACCGGG---GTC')


def test_put_gaps_in_nuc_seq2():
    seq = 'TTTAAAACCGGGGTC'
    trl = '-FKTG-V'

    eq_(sequence.put_gaps_in_nuc_seq(seq, trl), '---TTTAAAACCGGG---GTC')


def test_put_gaps_in_nuc_seq3():
    seq = 'TTTAAAACCGGGGTC'
    trl = 'FKTG-V-'

    eq_(sequence.put_gaps_in_nuc_seq(seq, trl), 'TTTAAAACCGGG---GTC---')


def test_put_gaps_in_nuc_seq4():
    seq = 'TTTAAAACCGGGGTC'
    trl = '-FKTG-V-'

    eq_(sequence.put_gaps_in_nuc_seq(seq, trl), '---TTTAAAACCGGG---GTC---')


def test_put_gaps_in_nuc_seq5():
    #trim sequence by default
    seq = 'TTTAAAACCGGGGTCAA'
    trl = '-FKTG-V-'

    eq_(sequence.put_gaps_in_nuc_seq(seq, trl), '---TTTAAAACCGGG---GTC---')


def test_put_gaps_in_nuc_seq6():
    #don't trim sequence by default
    seq = 'TTTAAAACCGGGGTCAA'
    trl = '-FKTG-V-'

    eq_(
        sequence.put_gaps_in_nuc_seq(seq, trl, False),
        '---TTTAAAACCGGG---GTC---AA'
    )


def test_get_seq_expected_syn_count():
    seq = 'ATGCATCGACTCTGCACTACG'

    eq_(sequence.get_seq_expected_syn_count(seq), (15, 48))


def test_alignment1():
    alg = sequence.Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCACT')
        ]
    )
    eq_(alg.get_position(6), 'TCA')


def test_alignment2():
    alg = sequence.Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCACT')
        ]
    )
    eq_(alg.get_seq_len(), 9)


def test_alignment3():
    alg = sequence.Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCCCT')
        ]
    )
    eq_(alg.get_consensus(), 'ACCCACCCT')


def test_alignment4():
    alg = sequence.Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCCCT')
        ]
    )
    ref_seq = alg.get_consensus()
    snps = [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')]
    eq_(alg.get_snps(ref_seq), snps)


def test_alignment5():
    alg = sequence.Alignment(
        seqs=[
            (1, 'ACTCACTCA'),
            (2, 'ACCCACCCT'),
            (3, 'ACCCGCCCT')
        ]
    )
    ref_seq = alg.get_consensus()
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
    eq_(alg.get_snps(ref_seq, full_size=True), snps)


def test_get_seq_number_of_syn1():
    ref_seq = 'ACCCACCCT'
    snps = [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')]
    counts = sequence.get_seq_number_of_syn(ref_seq, snps, start=0,
                                            trans_table=trans_tables.UNIVERSAL)
    eq_(counts, (2, 2))


def test_get_seq_number_of_syn2():
    ref_seq = 'ACCCACCCT'
    snps = [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')]
    counts = sequence.get_seq_number_of_syn(ref_seq, snps, start=1,
                                            trans_table=trans_tables.UNIVERSAL)
    eq_(counts, (1, 2))


def test_get_seq_number_of_syn3():
    ref_seq = 'ACCCACCCT'
    snps = [(2, 'T'), (4, 'G'), (6, 'T'), (8, 'A')]
    counts = sequence.get_seq_number_of_syn(ref_seq, snps, start=2,
                                            trans_table=trans_tables.UNIVERSAL)
    eq_(counts, (1, 2))
