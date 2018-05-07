import pytest
import random
import functools
from conftest import skip_no_connection

from mgkit.filter.taxon import filter_taxon_by_id_list
from mgkit.taxon import is_ancestor


pytestmark = skip_no_connection

def test_taxon_id_list1(ncbi_taxonomy):
    # at least prevotella as genus should be amog the others
    taxon_id = ncbi_taxonomy.find_by_name('prevotella')[0]
    filter_list = [
        taxon.taxon_id
        for taxon in ncbi_taxonomy
        if 'prevotella' in taxon.s_name.lower()
    ]
    assert filter_taxon_by_id_list(taxon_id, filter_list)


def test_taxon_id_list1_anc(ncbi_taxonomy):
    # at least prevotella as genus should be among the others
    taxon_id = random.choice(
        [
            taxon.taxon_id
            for taxon in ncbi_taxonomy
            if taxon.s_name.lower().startswith('prevotella ')
        ]
    )
    filter_list = ncbi_taxonomy.find_by_name('prevotella')
    assert filter_taxon_by_id_list(
        taxon_id,
        filter_list,
        func=functools.partial(is_ancestor, ncbi_taxonomy)
    )


def test_taxon_id_list2(ncbi_taxonomy):
    # no prevotella species will be found
    filter_list = ncbi_taxonomy.find_by_name('prevotella')
    taxon_id = random.choice(
        [
            taxon.taxon_id
            for taxon in ncbi_taxonomy
            if ('prevotella' in taxon.s_name) and (taxon.rank != 'genus')
        ]
    )
    assert not filter_taxon_by_id_list(taxon_id, filter_list)


def test_taxon_id_list2_anc(ncbi_taxonomy):
    # at least prevotella as genus should be among the others
    taxon_id = 839
    filter_list = ncbi_taxonomy.find_by_name('prevotella')
    assert not filter_taxon_by_id_list(
        taxon_id,
        filter_list,
        exclude=True,
        func=functools.partial(is_ancestor, ncbi_taxonomy)
    )
