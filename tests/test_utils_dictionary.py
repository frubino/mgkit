import pytest
from builtins import range
from mgkit.utils.dictionary import find_id_in_dict, combine_dict, \
    reverse_mapping, link_ids, combine_dict_one_value, merge_dictionaries, \
    split_dictionary_by_value, apply_func_to_values


@pytest.fixture
def valuedict():
    return {
        1: [1.0, 2.0],
        2: [2.0, 3.0],
        3: [4.0, 5.0, 6.0],
        4: [7.0]
    }


def test_combine_dict(valuedict):
    keydict = {
        'A': [1, 2, 3],
        'B': [3, 4]
    }
    result_dict = {
        'A': {1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
        'B': {4.0, 5.0, 6.0, 7.0},
    }
    assert combine_dict(keydict, valuedict) == result_dict


def test_combine_dict_one_value(valuedict):
    keydict = {
        'A': 1,
        'B': 2
    }
    result_dict = {
        'A': [1.0, 2.0],
        'B': [2.0, 3.0],
    }
    assert combine_dict_one_value(keydict, valuedict) == result_dict


@pytest.fixture
def id_map():
    return {
        'k1': ['v2', 'v3'],
        'k2': ['v3'],
        'k3': ['v4', 'v1'],
        'k4': ['v2', 'v1'],
    }


def test_link_ids1(id_map):

    result_dict = {
        'k3': set(['k4']),
        'k2': set(['k1']),
        'k1': set(['k2', 'k4']),
        'k4': set(['k3', 'k1'])
    }

    assert link_ids(id_map) == result_dict


# Try to parametrise this another time
def test_link_ids2(id_map):

    result_dict = {
        'k3': set(['k4']),
        'k2': set([]),
        'k1': set(['k4']),
        'k4': set(['k3', 'k1'])
    }

    assert link_ids(id_map, black_list=['v3']) == result_dict


def test_reverse_mapping(id_map):
    result_dict = {
        'v1': set(['k3', 'k4']),
        'v2': set(['k1', 'k4']),
        'v3': set(['k2', 'k1']),
        'v4': set(['k3'])
    }
    assert reverse_mapping(id_map) == result_dict


def test_find_id_in_dict(id_map):
    assert sorted(find_id_in_dict('v1', id_map)) == ['k3', 'k4']


def test_merge_dictionaries():

    res = merge_dictionaries(
        [
            {
                1: range(5)
            },
            {
                1: range(8)
            },
            {
                1: 9
            }
        ]
    )
    assert res == {1: {0, 1, 2, 3, 4, 5, 6, 7, 9}}


@pytest.mark.parametrize(
    "threshold,key_filter,result_dict1,result_dict2",
    [
        (6, None, {1: [1.0, 2.0], 2: [2.0, 3.0]}, {3: [4.0, 5.0, 6.0], 4: [7.0]}),
        (5, None, {1: [1.0, 2.0]}, {2: [2.0, 3.0], 3: [4.0, 5.0, 6.0], 4: [7.0]}),
        (5, [1, 2, 4], {1: [1.0, 2.0]}, {2: [2.0, 3.0], 4: [7.0]}),
    ]
)
def test_split_dictionary_by_value(threshold, key_filter, result_dict1, result_dict2):
    assert split_dictionary_by_value(valuedict(), threshold, sum, key_filter) == (result_dict1, result_dict2)


def test_apply_func_to_values():
    d = {
        1: ('vg', 'vt')
    }

    assert apply_func_to_values(d, lambda x: x[0]) == {1: {'v'}}
