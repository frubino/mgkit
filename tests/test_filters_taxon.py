from nose.tools import *
import random
import functools

from mgkit.filter.taxon import *
from mgkit.taxon import is_ancestor
# from _utils import skip_test

import taxon_data


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_lineage1():
    taxa = filter_taxonomy_by_lineage(taxon_data.TAXONOMY, 'archaea')
    ok_(
        all(
            'archaea' in taxon.lineage
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_lineage2():
    taxa = filter_taxonomy_by_lineage(taxon_data.TAXONOMY, 'bacteria')
    ok_(
        all(
            'bacteria' in taxon.lineage
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_lineage3():
    taxa = filter_taxonomy_by_lineage(taxon_data.TAXONOMY, 'archaea')
    ok_(
        all(
            'bacteria' not in taxon.lineage
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_lineage4():
    taxa = filter_taxonomy_by_lineage(taxon_data.TAXONOMY, 'bacteria')
    ok_(
        all(
            'archaea' not in taxon.lineage
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_rank1():
    taxa = filter_taxonomy_by_rank(taxon_data.TAXONOMY, 'genus')
    ok_(
        all(
            taxon.rank == 'genus'
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_rank2():
    taxa = filter_taxonomy_by_rank(taxon_data.TAXONOMY, 'genus')
    ok_(
        all(
            taxon.rank != 'order'
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_rank3():
    taxa = filter_taxonomy_by_rank(taxon_data.TAXONOMY, 'phylum')
    ok_(
        all(
            taxon.rank == 'phylum'
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_rank4():
    taxa = filter_taxonomy_by_rank(taxon_data.TAXONOMY, 'phylum')
    ok_(
        all(
            taxon.rank != 'order'
            for taxon in taxa
        )
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_taxon_id_list1():
    # at least prevotella as genus should be amog the others
    taxon_id = taxon_data.TAXONOMY.find_by_name('prevotella')[0]
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'prevotella' in taxon.s_name
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
            if taxon.s_name.startswith('prevotella ')
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
