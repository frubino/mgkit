import pytest
from conftest import skip_no_connection
from mgkit.snps.classes import GeneSNP
from mgkit.snps.filter import filter_genesyn_by_coverage, \
    filter_genesyn_by_gene_id, filter_genesyn_by_taxon_id, FilterFails
import mgkit.snps.mapper
import mgkit.taxon
import functools


@skip_no_connection
def test_snps_taxon1(ncbi_taxonomy):
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in ncbi_taxonomy
        if 'methanobrevibacter' in taxon.s_name.lower()
    ]
    taxon_id = ncbi_taxonomy.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSNP(taxon_id=taxon_id)


    assert filter_genesyn_by_taxon_id(
        gene_syn,
        ncbi_taxonomy,
        filter_list=filter_list,
        exclude=False
    )


@skip_no_connection
def test_snps_taxon1_rev(ncbi_taxonomy):
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in ncbi_taxonomy
        if 'methanobrevibacter' in taxon.s_name.lower()
    ]
    taxon_id = ncbi_taxonomy.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSNP(taxon_id=taxon_id)

    assert not filter_genesyn_by_taxon_id(
        gene_syn,
        ncbi_taxonomy,
        filter_list=filter_list,
        exclude=True
    )



@skip_no_connection
def test_snps_taxon2(ncbi_taxonomy):
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in ncbi_taxonomy
        if 'clostridium' in taxon.s_name
    ]
    taxon_id = ncbi_taxonomy.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSNP(taxon_id=taxon_id)

    assert not filter_genesyn_by_taxon_id(
        gene_syn,
        ncbi_taxonomy,
        filter_list=filter_list,
        exclude=False
    )


@skip_no_connection
def test_snps_taxon2_rev(ncbi_taxonomy):
    # will find it
    filter_list = [
        taxon.taxon_id
        for taxon in ncbi_taxonomy
        if 'clostridium' in taxon.s_name
    ]
    taxon_id = ncbi_taxonomy.find_by_name('methanobrevibacter')[0]
    gene_syn = GeneSNP(taxon_id=taxon_id)

    assert filter_genesyn_by_taxon_id(
        gene_syn,
        ncbi_taxonomy,
        filter_list=filter_list,
        exclude=True
    )


def test_snps_taxon3_exc1():
    gene_syn = GeneSNP()
    with pytest.raises(FilterFails):
        filter_genesyn_by_taxon_id(
            gene_syn, filter_list=range(10), taxonomy=None, func=filter
        )


def test_snps_taxon3_exc2():
    gene_syn = GeneSNP()
    with pytest.raises(FilterFails):
        filter_genesyn_by_taxon_id(
            gene_syn, filter_list=None
        )


def test_snps_gene_id1():
    gene_syn = GeneSNP(gene_id='K01201')
    gene_list = ['K01201', 'K02201', 'K01251']

    assert filter_genesyn_by_gene_id(
        gene_syn,
        gene_list,
        id_func=lambda x: x.gene_id
    )


def test_snps_gene_id2():
    gene_syn = GeneSNP(gene_id='K01201')
    gene_list = ['K01201', 'K02201', 'K01251']

    assert not filter_genesyn_by_gene_id(
        gene_syn,
        gene_list,
        exclude=True,
        id_func=lambda x: x.gene_id
    )


def test_snps_gene_id_exc():
    gene_syn = GeneSNP(gene_id='K01201')
    gene_list = None
    with pytest.raises(FilterFails):
        filter_genesyn_by_gene_id(gene_syn, gene_ids=gene_list)


def test_snps_gene_coverage1():
    gene_syn = GeneSNP(gene_id='K01201', coverage=4)
    min_cov = 4

    assert filter_genesyn_by_coverage(gene_syn, min_cov=min_cov)


def test_snps_gene_coverage2():
    gene_syn = GeneSNP(gene_id='K01201', coverage=3)
    min_cov = 4

    assert not filter_genesyn_by_coverage(gene_syn, min_cov=min_cov)


def test_snps_gene_coverage_exc():
    gene_syn = GeneSNP(gene_id='K01201')
    min_cov = None
    with pytest.raises(FilterFails):
        filter_genesyn_by_coverage(gene_syn, min_cov=min_cov)


def test_map_gene_id1():
    gene_map = {
        'K1': ['K2', 'K3']
    }
    assert list(mgkit.snps.mapper.map_gene_id('K1', gene_map)) == gene_map['K1']


def test_map_gene_id2():
    gene_map = {
        'K1': ['K2', 'K3']
    }

    assert list(mgkit.snps.mapper.map_gene_id('K2', gene_map)) == []


@skip_no_connection
def test_map_taxon_id_to_rank1(ncbi_taxonomy):
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_rank(
            839,
            rank='genus',
            taxonomy=ncbi_taxonomy,
            include_higher=False
        )
    ) == [838]


@skip_no_connection
def test_map_taxon_id_to_rank2(ncbi_taxonomy):
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_rank(
            2,
            rank='genus',
            taxonomy=ncbi_taxonomy,
            include_higher=False
        )
    ) == []


@skip_no_connection
def test_map_taxon_id_to_rank3(ncbi_taxonomy):
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_rank(
            2,
            rank='genus',
            taxonomy=ncbi_taxonomy,
            include_higher=True
        )
    ) == [2]


@skip_no_connection
def test_map_taxon_id_to_ancestor1(ncbi_taxonomy):
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        ncbi_taxonomy
    )
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_ancestor(
            839,
            anc_ids=[838],
            func=func
        )
    ) == [838]


@skip_no_connection
def test_map_taxon_id_to_ancestor2(ncbi_taxonomy):
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        ncbi_taxonomy
    )
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_ancestor(
            839,
            anc_ids=[838, 2],
            func=func
        )
    ) == [838, 2]


@skip_no_connection
def test_map_taxon_id_to_ancestor3(ncbi_taxonomy):
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        ncbi_taxonomy
    )
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_ancestor(
            839,
            anc_ids=[838, 2, 1485],
            func=func
        )
    ) == [838, 2]


@skip_no_connection
def test_map_taxon_id_to_ancestor4(ncbi_taxonomy):
    func = functools.partial(
        mgkit.taxon.is_ancestor,
        ncbi_taxonomy
    )
    assert list(
        mgkit.snps.mapper.map_taxon_id_to_ancestor(
            839,
            anc_ids=[2172, 1485],
            func=func
        )
    ) == []
