import pytest

from mgkit.snps.classes import GeneSNP, SNPType
from mgkit.snps.funcs import combine_sample_snps, flat_sample_snps

import numpy


def test_genesyn_calc_ratio1():
    # syn and nonsyn > 0
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    gs.add_snp(1, 'A', SNPType.syn)
    gs.add_snp(1, 'A', SNPType.syn)
    gs.add_snp(1, 'A', SNPType.nonsyn)
    gs.add_snp(1, 'A', SNPType.nonsyn)
    gs.add_snp(1, 'A', SNPType.nonsyn)
    assert gs.calc_ratio() == 1.0


def test_genesyn_calc_ratio2():
    # syn = nonsyn = 0, haplotypes=True
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    assert gs.calc_ratio(haplotypes=True) == 0.0


def test_genesyn_add1():
    # checks if the conserved gene_id and taxon_id is from the first instance
    gs1 = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6
    )
    gs2 = GeneSNP(
        exp_syn=4,
        exp_nonsyn=6,
        coverage=4
    )
    gs1.add(gs2)
    assert (gs1.gene_id, gs1.taxon_id) == ('A', 1)


def test_genesyn_add2():
    # checks if the values (syn, nonsyn, exp_syn, exp_nonsyn) are added
    gs1 = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6
    )
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.nonsyn)
    gs2 = GeneSNP(
        exp_syn=3,
        exp_nonsyn=5,
        coverage=4
    )
    gs2.add_snp(1, 'A', SNPType.syn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)

    gs1.add(gs2)
    assert (gs1.exp_syn, gs1.exp_nonsyn, gs1.syn, gs1.nonsyn) == (7, 11, 4, 3)


def test_genesyn_add3():
    # checks if the values if coverage is correctly added
    # coverage in gs2 is not None, in gs1 is None
    gs1 = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6
    )
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.nonsyn)
    gs2 = GeneSNP(
        exp_syn=3,
        exp_nonsyn=5,
        coverage=5
    )
    gs2.add_snp(1, 'A', SNPType.syn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)

    gs1.add(gs2)
    assert gs1.coverage == 5


def test_genesyn_add4():
    # checks if the values if coverage is correctly added
    # coverage in gs2 is not None, in gs1 is not None
    gs1 = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
        coverage=2
    )
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.nonsyn)
    gs2 = GeneSNP(
        exp_syn=3,
        exp_nonsyn=5,
        coverage=5
    )
    gs2.add_snp(1, 'A', SNPType.syn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)

    gs1.add(gs2)
    assert gs1.coverage == 7


def test_genesyn_add5():
    # checks if the values of coverage is correctly added
    # coverage in gs2 is None, in gs1 is None
    gs1 = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.syn)
    gs1.add_snp(1, 'A', SNPType.nonsyn)
    gs2 = GeneSNP(
        exp_syn=3,
        exp_nonsyn=5,
        coverage=None
    )
    gs2.add_snp(1, 'A', SNPType.syn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)
    gs2.add_snp(1, 'A', SNPType.nonsyn)

    gs1.add(gs2)
    assert gs1.coverage is None


def test_genesyn_calc_ratio3b():
    # syn = nonsyn = 0
    # coverage is None
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    assert numpy.isnan(gs.calc_ratio())


def test_genesyn_calc_ratio3c():
    # syn = nonsyn = 0
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    assert numpy.isnan(gs.calc_ratio())


def test_genesyn_calc_ratio3d():
    # syn = nonsyn = 0, haplotypes = True
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    assert gs.calc_ratio(haplotypes=True) == 0


def test_genesyn_calc_ratio4a():
    # nonsyn > 0, syn = 0
    # flag_value=True
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    gs.add_snp(1, 'A', SNPType.nonsyn)
    assert gs.calc_ratio_flag() == -1


def test_genesyn_calc_ratio4b():
    # nonsyn = 0, syn > 0
    # flag_value=True
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    gs.add_snp(1, 'A', SNPType.syn)
    assert gs.calc_ratio_flag() == -2


def test_genesyn_calc_ratio4c():
    # nonsyn = 0, syn > 0
    # flag_value=True
    gs = GeneSNP(
        gene_id='A',
        taxon_id=1,
        exp_syn=4,
        exp_nonsyn=6,
    )
    assert gs.calc_ratio_flag() == -3


SNP_DATA = {
    'sample1': {
        'gene1': GeneSNP(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=6,
            exp_nonsyn=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        ),  # pN/pS = 1.0
        'gene2': GeneSNP(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        )  # pN/pS = 2.0
    },
    'sample2': {
        'gene1': GeneSNP(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=3,
            exp_nonsyn=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        ),  # pN/pS = 1.0
        'gene2': GeneSNP(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=6,
            exp_nonsyn=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        ),  # pN/pS = 2.0
        'gene3': GeneSNP(
            gene_id='gene3',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        )  # pN/pS = 2.0
    }

}


SNP_DATA2 = {
    'sample1': {
        'gene1': GeneSNP(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=6,
            exp_nonsyn=4,
            coverage=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        ),  # pN/pS = 1.0
        'gene2': GeneSNP(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            coverage=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        )  # pN/pS = 2.0
    },
    'sample2': {
        'gene1': GeneSNP(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=3,
            exp_nonsyn=4,
            coverage=4,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        ),  # pN/pS = 1.0
        'gene2': GeneSNP(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=6,
            exp_nonsyn=4,
            coverage=3,
            snps=[
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.syn),
                (1, 'A', SNPType.nonsyn),
                (1, 'A', SNPType.nonsyn),
            ]
        )  # pN/pS = 2.0
    }

}


def test_flat_sample_snps1():
    data = flat_sample_snps(SNP_DATA2, 4)
    assert (
            data['all_samples']['gene1'].syn,
            data['all_samples']['gene2'].syn
        ) == (7, 3)


def test_flat_sample_snps2():
    data = flat_sample_snps(SNP_DATA2, 3)
    assert (
            data['all_samples']['gene1'].syn,
            data['all_samples']['gene2'].syn
        ) == (7, 8)


def test_combine_sample_min_num1():
    # check if min_num is working
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    df.shape == (3, 2)


def test_combine_sample_snps_min_num2():
    # check if min_num is working
    df = combine_sample_snps(SNP_DATA, 2, [lambda x: x])
    assert df.shape == (2, 2)


def test_combine_sample_snps_values():
    # check if values are correct is working
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    assert (df.min().min(), df.max().max()) == (0.5, 1.)


def test_combine_sample_snps_index1():
    # check index_type
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    assert df.index.tolist() == [('gene1', 839), ('gene2', 838), ('gene3', 838)]


def test_combine_sample_snps_index2():
    # check index_type
    df = combine_sample_snps(
        SNP_DATA,
        1,
        [lambda x: x],
        index_type='gene'
    )
    # Can't expect the list to be sorted
    assert sorted(df.index.tolist()) == ['gene1', 'gene2', 'gene3']


def test_combine_sample_snps_index3():
    # check index_type
    df = combine_sample_snps(
        SNP_DATA,
        1,
        [lambda x: x],
        index_type='taxon'
    )
    assert df.index.tolist() == [838, 839]
