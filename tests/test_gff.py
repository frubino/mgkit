from nose.tools import eq_, ok_, with_setup, raises

from mgkit.io import gff
from mgkit.utils import sequence
import misc_data


def test_genomicrange_contains_annotation1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.Annotation(start=25, end=30) in gen_range1,
        True
    )


def test_genomicrange_contains_annotation2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.Annotation(start=30, end=25) in gen_range1,
        True
    )


def test_genomicrange_contains_annotation3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.Annotation(start=19, end=30) in gen_range1,
        False
    )


def test_genomicrange_contains_annotation4():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.Annotation(start=25, end=31) in gen_range1,
        False
    )


def test_genomicrange_get_relative_pos1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(20),
        1
    )


def test_genomicrange_get_relative_pos2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(30),
        11
    )


def test_genomicrange_get_relative_pos3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(25),
        6
    )


@raises(ValueError)
def test_genomicrange_get_relative_pos_fail1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    gen_range1.get_relative_pos(19)


@raises(ValueError)
def test_genomicrange_get_relative_pos_fail2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    gen_range1.get_relative_pos(31)


def test_genomicrange_union1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 20)
    )


def test_genomicrange_union3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 20
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union4():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 21
    gen_range2.end = 30

    gen_range_u = gen_range2.union(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union_fail1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq2'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_union_fail2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '-'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_union_fail3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '-'
    gen_range2.start = 21
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (19, 20)
    )


def test_genomicrange_intersect2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 15
    gen_range2.end = 30

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (15, 20)
    )


def test_genomicrange_intersect3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 15
    gen_range2.end = 30

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (15, 20)
    )


def test_genomicrange_intersect4():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 20)
    )


def test_genomicrange_intersect5():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        len(gen_range_u),
        len(gen_range2)
    )


def test_genomicrange_intersect6():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        len(gen_range_u),
        len(gen_range2)
    )


def test_genomicrange_intersect_fail1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq2'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect_fail2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '-'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect_fail3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 30
    gen_range2.end = 40

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def elongate_data():
    seq_id = 'test1'

    test_ann = [
        gff.GenomicRange(seq_id=seq_id, start=1, end=10),
        gff.GenomicRange(seq_id=seq_id, start=10, end=15),
        gff.GenomicRange(seq_id=seq_id, start=16, end=18),
    ]
    return test_ann


def test_annotation1():
    ann = gff.Annotation(start=1, end=10)
    eq_(len(ann), 10)


def test_genomicrange_misc1():
    ann = gff.GenomicRange(start=1, end=10)
    eq_(len(ann), 10)
