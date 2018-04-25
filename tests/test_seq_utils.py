from __future__ import division
from nose.tools import eq_, ok_

import numpy
from mgkit.utils import sequence
from mgkit.utils import trans_tables


def test_sequence_composition3():
    seq = 'A' * 10 + 'C' * 4 + 'T' * 5 + 'G' * 11 + 'N' * 2
    eq_(
        sorted(sequence.sequence_composition(seq, chars=('A', 'C'))),
        [('A', 10), ('C', 4)]
    )


def test_get_kmers1():
    seq = 'ACTG' * 2
    eq_(
        list(sequence._get_kmers(seq, 4)),
        ['ACTG', 'CTGA', 'TGAC', 'GACT', 'ACTG']
    )


def test_get_kmers2():
    seq = 'ACTG' * 2
    eq_(
        list(sequence._get_kmers(seq, 5)),
        ['ACTGA', 'CTGAC', 'TGACT', 'GACTG']
    )


def test_sliding_window1():
    seq = 'ACTG' * 5
    eq_(
        list(
             sequence._sliding_window(seq, 4, 4)
        ),
        ['ACTG'] * 5
    )


def test_sliding_window2():
    seq = 'ACTG' * 2
    eq_(
        list(
             sequence._sliding_window(seq, 4, 2)
        ),
        ['ACTG', 'TGAC', 'ACTG']
    )


def test_sliding_window3():
    seq = 'ACTG' * 2
    eq_(
        list(
             sequence._sliding_window(seq, 4, 3)
        ),
        ['ACTG', 'GACT']
    )


def test_sequence_signature1():
    seq = 'ACTG' * 2
    count = sequence._sequence_signature(
        seq,
        4,
        4,
        4
    )
    eq_(
        len(count),
        2
    )


def test_sequence_signature2():
    seq = 'ACTG' * 2
    count = sequence._sequence_signature(
        seq,
        4,
        4,
        4
    )
    eq_(
        count[0]['ACTG'],
        1
    )


def test_sequence_signature3():
    seq = 'ACTG' * 2
    count = sequence._sequence_signature(
        seq,
        8,
        4,
        4
    )
    eq_(
        count[0]['ACTG'],
        2
    )


def test_sequence_signature4():
    seq = 'ACTG' * 2
    count = sequence._sequence_signature(
        seq,
        5,
        4,
        5
    )
    eq_(
        len(count),
        1
    )


def test_sequence_signature_cython():
    seq = 'ACTG' * 2
    countP = sequence._sequence_signature(
        seq,
        4,
        4,
        4
    )
    countC = sequence.sequence_signature(
        seq,
        4,
        4,
        4
    )
    eq_(
        len(countP),
        len(countC),
    )


def test_sliding_window_cython():
    seq = 'ACTG' * 5
    eq_(
        list(
             sequence._sliding_window(seq, 4, 4)
        ),
        list(
             sequence.sliding_window(seq, 4, 4)
        )
    )


def test_get_kmers_cython():
    seq = 'ACTG' * 2
    eq_(
        list(sequence._get_kmers(seq, 4)),
        list(sequence.get_kmers(seq, 4))
    )
