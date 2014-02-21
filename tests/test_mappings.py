from nose.tools import *

import functools
from mgkit.mappings.taxon import *
from mgkit.taxon import is_ancestor

import taxon_data


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_by_id_list1():
    result = map_taxon_by_id_list(
        839,
        [838, 1485],
        functools.partial(
            is_ancestor,
            taxon_data.TAXONOMY
        )
    )
    eq_(list(result), [838])


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_by_id_list2():
    result = map_taxon_by_id_list(
        2172,
        [838, 1485],
        functools.partial(
            is_ancestor,
            taxon_data.TAXONOMY
        )
    )
    eq_(list(result), [])
