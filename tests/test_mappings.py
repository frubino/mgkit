import pytest
from conftest import skip_no_connection, ncbi_taxonomy

import functools
from mgkit.mappings.taxon import map_taxon_by_id_list
from mgkit.taxon import is_ancestor


@skip_no_connection
def test_map_taxon_by_id_list1(ncbi_taxonomy):
    result = map_taxon_by_id_list(
        839,
        [838, 1485],
        functools.partial(
            is_ancestor,
            ncbi_taxonomy
        )
    )
    assert list(result) == [838]


@skip_no_connection
def test_map_taxon_by_id_list2(ncbi_taxonomy):
    result = map_taxon_by_id_list(
        2172,
        [838, 1485],
        functools.partial(
            is_ancestor,
            ncbi_taxonomy
        )
    )
    assert list(result) == []
