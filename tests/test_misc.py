from nose.tools import *

from mgkit.utils.common import average_length, between, union_range, \
    range_intersect
from mgkit.utils.dictionary import *


def test_avg_len1():
    eq_(average_length(1, 200, 3, 132), 165.0)


def test_avg_len2():
    eq_(average_length(4, 2245, 3, 13223), 7731.5)


def test_between1():
    eq_(between(1, 0, 10), True)


def test_between2():
    eq_(between(1, 2, 10), False)


def test_between3():
    eq_(between(1, 1, 10), True)


def test_between4():
    eq_(between(11, 1, 10), False)


def test_between5():
    eq_(between(10, 1, 10), True)


def test_combine_dict1():
    keydict = {
        'A': [1, 2, 3],
        'B': [3, 4]
    }
    valuedict = {
        1: [1.0, 2.0],
        2: [2.0, 3.0],
        3: [4.0, 5.0, 6.0],
        4: [7.0]
    }
    eq_(
        combine_dict(keydict, valuedict),
        {
            'A': {1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
            'B': {4.0, 5.0, 6.0, 7.0},
        }
    )


def test_combine_dict_one_value():
    keydict = {
        'A': 1,
        'B': 2
    }
    valuedict = {
        1: [1.0, 2.0],
        2: [2.0, 3.0],
        3: [4.0, 5.0, 6.0],
        4: [7.0]
    }
    eq_(
        combine_dict_one_value(keydict, valuedict),
        {
            'A': [1.0, 2.0],
            'B': [2.0, 3.0],
        }
    )


def test_link_ids1():
    id_map = {
        'k1': ['v2', 'v3'],
        'k2': ['v3'],
        'k3': ['v4', 'v1'],
        'k4': ['v2', 'v1'],
    }
    eq_(
        link_ids(id_map),
        {
            'k3': set(['k4']),
            'k2': set(['k1']),
            'k1': set(['k2', 'k4']),
            'k4': set(['k3', 'k1'])
        }
    )


def test_link_ids2():
    id_map = {
        'k1': ['v2', 'v3'],
        'k2': ['v3'],
        'k3': ['v4', 'v1'],
        'k4': ['v2', 'v1'],
    }
    eq_(
        link_ids(id_map, black_list=['v3']),
        {
            'k3': set(['k4']),
            'k2': set([]),
            'k1': set(['k4']),
            'k4': set(['k3', 'k1'])
        }
    )


def test_reverse_mapping():
    map_dict = {
        'k1': ['v2', 'v3'],
        'k2': ['v3'],
        'k3': ['v4', 'v1'],
        'k4': ['v2', 'v1'],
    }
    eq_(
        reverse_mapping(map_dict),
        {
            'v1': set(['k3', 'k4']),
            'v2': set(['k1', 'k4']),
            'v3': set(['k2', 'k1']),
            'v4': set(['k3'])
        }
    )


def test_find_id_in_dict():
    map_dict = {
        'k1': ['v2', 'v3'],
        'k2': ['v3'],
        'k3': ['v4', 'v1'],
        'k4': ['v2', 'v1'],
    }
    eq_(
        sorted(find_id_in_dict('v1', map_dict)),
        ['k3', 'k4']
    )


def test_union_range1():
    eq_(
        union_range(1, 10, 2, 13),
        (1, 13)
    )


def test_union_range2():
    eq_(
        union_range(1, 10, 10, 13),
        (1, 13)
    )


def test_union_range3():
    eq_(
        union_range(10, 13, 1, 10),
        (1, 13)
    )


def test_union_range4():
    eq_(
        union_range(1, 10, 11, 13),
        None
    )


def test_union_range5():
    eq_(
        union_range(10.0, 13.0, 1, 10),
        (1, 13)
    )


def test_range_intersect1():
    range1 = (10, 20)
    range2 = (19, 30)
    eq_(
        range_intersect(*(range1 + range2)),
        (19, 20)
    )


def test_range_intersect2():
    range1 = (10, 20)
    range2 = (15, 30)
    eq_(
        range_intersect(*(range1 + range2)),
        (15, 20)
    )


def test_range_intersect3():
    range1 = (10, 20)
    range2 = (10, 20)
    eq_(
        range_intersect(*(range1 + range2)),
        (10, 20)
    )


def test_range_intersect4():
    range1 = (10, 20)
    range2 = (12, 18)
    eq_(
        range_intersect(*(range1 + range2)),
        (12, 18)
    )


def test_range_intersect_fail1():
    range1 = (10, 20)
    range2 = (30, 40)
    eq_(
        range_intersect(*(range1 + range2)),
        None
    )
