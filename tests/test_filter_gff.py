import pytest
import functools

from mgkit.io.gff import Annotation
from mgkit.filter.gff import choose_annotation, filter_annotations, \
    filter_base, filter_len, filter_attr_num, filter_attr_str, filter_base_num


def test_choose_annotation_contained1():
    # same size, better score
    a = Annotation(start=1, end=10, bitscore=10)
    b = Annotation(start=1, end=10, bitscore=15)
    assert choose_annotation(a, b) == a


def test_choose_annotation_contained2():
    # same score, longer
    a = Annotation(start=1, end=10, bitscore=10)
    b = Annotation(start=5, end=10, bitscore=10)
    assert choose_annotation(a, b) == b


def test_choose_annotation_overlap1():
    # overlapping, better score
    a = Annotation(start=1, end=10, bitscore=10)
    b = Annotation(start=5, end=11, bitscore=15)
    assert choose_annotation(a, b, overlap=4) == a


def test_choose_annotation_overlap2():
    # overlapping, same score, choose longer
    a = Annotation(start=1, end=11, bitscore=10)
    b = Annotation(start=5, end=11, bitscore=10)
    assert choose_annotation(a, b, overlap=4) == b


def test_choose_annotation_none1():
    # no overlap
    a = Annotation(start=1, end=10, bitscore=10)
    b = Annotation(start=11, end=12, bitscore=10)
    assert choose_annotation(a, b, overlap=4) is None


def test_choose_annotation_none2():
    # overlap lower than the limits
    a = Annotation(start=1, end=10, bitscore=10)
    b = Annotation(start=9, end=12, bitscore=10)
    assert choose_annotation(a, b, overlap=2) is None


def test_filter_annotations1():

    annotations = [
        Annotation(start=1, end=10, bitscore=20),
        Annotation(start=5, end=10, bitscore=30),
        Annotation(start=9, end=15, bitscore=10),
        Annotation(start=16, end=20, bitscore=30),
        Annotation(start=16, end=20, bitscore=20),
        Annotation(start=21, end=30, bitscore=40),
    ]
    choose_func = functools.partial(choose_annotation, overlap=2)

    assert filter_annotations(annotations, choose_func=choose_func) == \
        set(
            [
                annotations[1],
                annotations[2],
                annotations[3],
                annotations[5],
            ]
        )


def test_filter_base_ok1():
    a = Annotation(seq_id='seq1')
    assert filter_base(a, 'seq_id', 'seq1')


def test_filter_base_fail1():
    a = Annotation(seq_id='seq1')
    assert not filter_base(a, 'seq_id', 'seq')


def test_filter_base_fail2():
    a = Annotation(seq_id='seq1')
    assert not filter_base(a, 'seqid', 'seq1')


def test_filter_len_ge_ok1():
    # equal
    a = Annotation(start=1, end=2)
    assert filter_len(a, 2, True)


def test_filter_len_le_ok1():
    # equal
    a = Annotation(start=1, end=2)
    assert filter_len(a, 2, False)


def test_filter_len_ge_ok2():
    # greater
    a = Annotation(start=1, end=3)
    assert filter_len(a, 2, True)


def test_filter_len_le_ok2():
    # lower
    a = Annotation(start=1, end=3)
    assert filter_len(a, 5, False)


def test_filter_len_ge_fail1():
    a = Annotation(start=1, end=3)
    assert not filter_len(a, 5, True)


def test_filter_len_le_fail1():
    a = Annotation(start=1, end=10)
    assert not filter_len(a, 5, False)


def test_filter_base_num_ge_ok1():
    a = Annotation(bitscore=10)
    assert filter_base_num(a, 'bitscore', 5, True)


def test_filter_base_num_ge_ok2():
    a = Annotation(bitscore=10)
    assert filter_base_num(a, 'bitscore', 10, True)


def test_filter_base_num_ge_fail1():
    a = Annotation(bitscore=10)
    assert not filter_base_num(a, 'bitscore', 15, True)


def test_filter_base_num_ge_fail2():
    a = Annotation()
    assert not filter_base_num(a, 'bitscore', 15, True)


def test_filter_base_num_le_ok1():
    a = Annotation(bitscore=10)
    assert filter_base_num(a, 'bitscore', 15, False)


def test_filter_base_num_le_ok2():
    a = Annotation(bitscore=10)
    assert filter_base_num(a, 'bitscore', 10, False)


def test_filter_base_num_le_fail1():
    a = Annotation(bitscore=10)
    assert not filter_base_num(a, 'bitscore', 5, False)


def test_filter_base_num_le_fail2():
    a = Annotation()
    assert not filter_attr_num(a, 'bitscore', 15, False)


def test_filter_attr_num_ge_ok1():
    a = Annotation(bitscore=10)
    assert filter_attr_num(a, 'bitscore', 5, True)


def test_filter_attr_num_ge_ok2():
    a = Annotation(bitscore=10)
    assert filter_attr_num(a, 'bitscore', 10, True)


def test_filter_attr_num_ge_fail1():
    a = Annotation(bitscore=10)
    assert not filter_attr_num(a, 'bitscore', 15, True)


def test_filter_attr_num_ge_fail2():
    a = Annotation()
    assert not filter_attr_num(a, 'bitscore', 15, True)


def test_filter_attr_num_le_ok1():
    a = Annotation(bitscore=10)
    assert filter_attr_num(a, 'bitscore', 15, False)


def test_filter_attr_num_le_ok2():
    a = Annotation(bitscore=10)
    assert filter_attr_num(a, 'bitscore', 10, False)


def test_filter_attr_num_le_fail1():
    a = Annotation(bitscore=10)
    assert not filter_attr_num(a, 'bitscore', 5, False)


def test_filter_attr_num_le_fail2():
    a = Annotation()
    assert not filter_attr_num(a, 'bitscore', 5, False)


def test_filter_attr_str_eq_ok1():
    a = Annotation(ec='1.1.1.1')
    assert filter_attr_str(a, 'ec', '1.1.1.1', True)


def test_filter_attr_str_eq_fail1():
    a = Annotation(ec='1.1.1.1')
    assert not filter_attr_str(a, 'ec1', '1.1.1.1', True)


def test_filter_attr_str_eq_fail2():
    a = Annotation(ec='1.1.1.1')
    assert not filter_attr_str(a, 'ec', '1.1.1', True)


def test_filter_attr_str_in_ok1():
    a = Annotation(ec='1.1.1.1')
    assert filter_attr_str(a, 'ec', '1.1.1.1', False)


def test_filter_attr_str_in_ok2():
    a = Annotation(ec='1.1.1.1')
    assert filter_attr_str(a, 'ec', '1.1.1', False)


def test_filter_attr_str_in_fail1():
    a = Annotation(ec='1.1.1.1')
    assert not filter_attr_str(a, 'ec1', '1.1.1.1', False)


def test_filter_attr_str_in_fail2():
    a = Annotation(ec='1.1.1.1')
    assert not filter_attr_str(a, 'ec', '1.1.1.2', False)
