from nose.tools import eq_

from mgkit.utils.dictionary import find_id_in_dict, combine_dict, \
    reverse_mapping, link_ids, combine_dict_one_value, merge_dictionaries


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


def test_merge_dictionaries():

    res = merge_dictionaries(
        [
            {
                1: xrange(5)
            },
            {
                1: range(8)
            },
            {
                1: 9
            }
        ]
    )
    eq_(
        res,
        {1: {0, 1, 2, 3, 4, 5, 6, 7, 9}}
    )
