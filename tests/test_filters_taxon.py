from nose.tools import ok_, eq_, with_setup
import random
import functools

from mgkit.filter.taxon import filter_taxon_by_id_list
from mgkit.taxon import is_ancestor
# from _utils import skip_test

import taxon_data

@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_id_list1():
    # at least prevotella as genus should be amog the others
    taxon_id = taxon_data.TAXONOMY.find_by_name('prevotella')[0]
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'prevotella' in taxon.s_name.lower()
    ]
    eq_(
        filter_taxon_by_id_list(taxon_id, filter_list), True
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_id_list1_anc():
    # at least prevotella as genus should be among the others
    taxon_id = random.choice(
        [
            taxon.taxon_id
            for taxon in taxon_data.TAXONOMY
            if taxon.s_name.lower().startswith('prevotella ')
        ]
    )
    filter_list = taxon_data.TAXONOMY.find_by_name('prevotella')
    eq_(
        filter_taxon_by_id_list(
            taxon_id,
            filter_list,
            func=functools.partial(is_ancestor, taxon_data.TAXONOMY)
        ),
        True
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_id_list2():
    # no prevotella species will be found
    filter_list = taxon_data.TAXONOMY.find_by_name('prevotella')
    taxon_id = random.choice(
        [
            taxon.taxon_id
            for taxon in taxon_data.TAXONOMY
            if ('prevotella' in taxon.s_name) and (taxon.rank != 'genus')
        ]
    )
    eq_(
        filter_taxon_by_id_list(taxon_id, filter_list), False
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_id_list2_anc():
    # at least prevotella as genus should be among the others
    taxon_id = 839
    filter_list = taxon_data.TAXONOMY.find_by_name('prevotella')
    eq_(
        filter_taxon_by_id_list(
            taxon_id,
            filter_list,
            exclude=True,
            func=functools.partial(is_ancestor, taxon_data.TAXONOMY)
        ),
        False
    )
