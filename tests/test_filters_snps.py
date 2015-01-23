from nose.tools import *

from mgkit.snps.classes import GeneSyn
from mgkit.snps.filter import *
import mgkit.snps.mapper
import mgkit.taxon
import functools

import taxon_data


@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon1():
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'methanobrevibacter' in taxon.s_name
    ]
    taxon_id = taxon_data.TAXONOMY.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSyn(taxon_id=taxon_id)

    eq_(
        filter_genesyn_by_taxon_id(
            gene_syn,
            taxon_data.TAXONOMY,
            filter_list=filter_list,
            exclude=False
        ),
        True
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon1_rev():
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'methanobrevibacter' in taxon.s_name
    ]
    taxon_id = taxon_data.TAXONOMY.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSyn(taxon_id=taxon_id)

    eq_(
        filter_genesyn_by_taxon_id(
            gene_syn,
            taxon_data.TAXONOMY,
            filter_list=filter_list,
            exclude=True
        ),
        False
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon2():
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'clostridium' in taxon.s_name
    ]
    taxon_id = taxon_data.TAXONOMY.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSyn(taxon_id=taxon_id)

    eq_(
        filter_genesyn_by_taxon_id(
            gene_syn,
            taxon_data.TAXONOMY,
            filter_list=filter_list,
            exclude=False
        ),
        False
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon2_rev():
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in taxon_data.TAXONOMY
        if 'clostridium' in taxon.s_name
    ]
    taxon_id = taxon_data.TAXONOMY.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSyn(taxon_id=taxon_id)

    eq_(
        filter_genesyn_by_taxon_id(
            gene_syn,
            taxon_data.TAXONOMY,
            filter_list=filter_list,
            exclude=True
        ),
        True
    )


@raises(FilterFails)
@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon3_exc1():
    gene_syn = GeneSyn()
    filter_genesyn_by_taxon_id(
        gene_syn, filter_list=range(10), taxonomy=None, func=filter
    )


@raises(FilterFails)
@with_setup(setup=taxon_data.setup_taxon_data)
def test_snps_taxon3_exc2():
    gene_syn = GeneSyn()
    filter_genesyn_by_taxon_id(
        gene_syn, filter_list=None
    )


def test_snps_gene_id1():
    gene_syn = GeneSyn(gene_id='K01201')
    gene_list = ['K01201', 'K02201', 'K01251']

    eq_(
        filter_genesyn_by_gene_id(
            gene_syn,
            gene_list,
            id_func=lambda x: x.gene_id
        ),
        True
    )


def test_snps_gene_id2():
    gene_syn = GeneSyn(gene_id='K01201')
    gene_list = ['K01201', 'K02201', 'K01251']

    eq_(
        filter_genesyn_by_gene_id(
            gene_syn,
            gene_list,
            exclude=True,
            id_func=lambda x: x.gene_id
        ),
        False
    )


@raises(FilterFails)
def test_snps_gene_id_exc():
    gene_syn = GeneSyn(gene_id='K01201')
    gene_list = None
    filter_genesyn_by_gene_id(gene_syn, gene_ids=gene_list)


def test_snps_gene_coverage1():
    gene_syn = GeneSyn(gene_id='K01201', coverage=4)
    min_cov = 4

    eq_(filter_genesyn_by_coverage(gene_syn, min_cov=min_cov), True)


def test_snps_gene_coverage2():
    gene_syn = GeneSyn(gene_id='K01201', coverage=3)
    min_cov = 4

    eq_(filter_genesyn_by_coverage(gene_syn, min_cov=min_cov), False)


@raises(FilterFails)
def test_snps_gene_coverage_exc():
    gene_syn = GeneSyn(gene_id='K01201')
    min_cov = None
    filter_genesyn_by_coverage(gene_syn, min_cov=min_cov)


def test_map_gene_id1():
    gene_map = {
        'K1': ['K2', 'K3']
    }
    eq_(
        list(mgkit.snps.mapper.map_gene_id('K1', gene_map)),
        gene_map['K1']
    )


def test_map_gene_id2():
    gene_map = {
        'K1': ['K2', 'K3']
    }
    eq_(
        list(mgkit.snps.mapper.map_gene_id('K2', gene_map)),
        []
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_rank1():
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_rank(
                839,
                rank='genus',
                taxonomy=taxon_data.TAXONOMY,
                include_higher=False
            )
        ),
        [838]
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_rank2():
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_rank(
                2,
                rank='genus',
                taxonomy=taxon_data.TAXONOMY,
                include_higher=False
            )
        ),
        []
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_rank3():
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_rank(
                2,
                rank='genus',
                taxonomy=taxon_data.TAXONOMY,
                include_higher=True
            )
        ),
        [2]
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_ancestor1():
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        taxon_data.TAXONOMY
    )
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_ancestor(
                839,
                anc_ids=[838],
                func=func
            )
        ),
        [838]
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_ancestor2():
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        taxon_data.TAXONOMY
    )
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_ancestor(
                839,
                anc_ids=[838, 2],
                func=func
            )
        ),
        [838, 2]
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_ancestor3():
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        taxon_data.TAXONOMY
    )
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_ancestor(
                839,
                anc_ids=[838, 2, 1485],
                func=func
            )
        ),
        [838, 2]
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_map_taxon_id_to_ancestor4():
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        taxon_data.TAXONOMY
    )
    eq_(
        list(
            mgkit.snps.mapper.map_taxon_id_to_ancestor(
                839,
                anc_ids=[2172, 1485],
                func=func
            )
        ),
        []
    )
