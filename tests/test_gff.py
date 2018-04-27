from nose.tools import eq_, ok_, with_setup, raises

from mgkit.io import gff
from mgkit.utils import sequence
import misc_data


@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_from_hmmer():
    checks = (
        'contig-1442648',
        'K00001_4479_poaceae',
        693,
        894,
        'K00001',
        4479,
        'poaceae'
    )
    ann = gff.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS
    )

    eq_(
        (
            ann.seq_id,
            ann.attr['name'],
            ann.attr['aa_from'],
            ann.attr['aa_to'],
            ann.gene_id,
            ann.taxon_id,
            ann.attr['taxon_name']
        ),
        checks
    )


def test_gff_glimmer3_line1():
    header = 'sequence0001'
    line = 'orf00001       66      611  +3     6.08'
    annotation = gff.from_glimmer3(header, line)

    eq_(
        (
            annotation.seq_id,
            annotation.start,
            annotation.end,
            annotation.score,
            annotation.strand,
            annotation.phase,
            annotation.attr['orf_id'],
            annotation.attr['frame']
        ),
        (
            header,
            66,
            611,
            6.08,
            '+',
            2,
            'orf00001',
            '+3'
        )
    )


def test_gff_glimmer3_line2():
    header = 'sequence0001'
    line = 'orf00001       66      11  -2     6.08'
    annotation = gff.from_glimmer3(header, line)

    eq_(
        (
            annotation.start,
            annotation.end,
            annotation.strand,
            annotation.phase,
            annotation.attr['frame']
        ),
        (
            11,
            66,
            '-',
            1,
            '-2'
        )
    )


@with_setup(setup=misc_data.setup_nucseq_data)
def test_gff_from_sequence1():
    annotation = gff.from_sequence(
        'contig-110637',
        misc_data.NUC_SEQS['contig-110637']
    )
    eq_(
        (
            annotation.seq_id,
            annotation.start,
            annotation.end,
        ),
        (
            'contig-110637',
            1,
            len(misc_data.NUC_SEQS['contig-110637'])
        )
    )


def test_genomicrange_contains1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        1 in gen_range1,
        False
    )


def test_genomicrange_contains2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        20 in gen_range1,
        True
    )


def test_genomicrange_contains3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        30 in gen_range1,
        True
    )


def test_genomicrange_contains4():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        31 in gen_range1,
        False
    )


def test_genomicrange_contains_tuple1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        (25, 30) in gen_range1,
        True
    )


def test_genomicrange_contains_tuple2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        (25, 31) in gen_range1,
        False
    )


def test_genomicrange_contains_tuple3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        (19, 30) in gen_range1,
        False
    )


def test_genomicrange_contains_tuple4():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        (30, 25) in gen_range1,
        True
    )


def test_genomicrange_contains_genomicrange1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.GenomicRange(start=25, end=30) in gen_range1,
        True
    )


def test_genomicrange_contains_genomicrange2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.GenomicRange(start=30, end=25) in gen_range1,
        True
    )


def test_genomicrange_contains_genomicrange3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.GenomicRange(start=19, end=30) in gen_range1,
        False
    )


def test_genomicrange_contains_genomicrange4():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gff.GenomicRange(start=25, end=31) in gen_range1,
        False
    )


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
