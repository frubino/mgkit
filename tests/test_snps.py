from nose.tools import *

from mgkit.snps.classes import GeneSyn
from mgkit.snps.funcs import combine_sample_snps, flat_sample_snps

import numpy


def test_genesyn_init1():
    gs = GeneSyn(
        gid='A',
        taxon=1
    )
    eq_(
       (gs.gid, gs.taxon),
       (gs.gene_id, gs.taxon_id)
    )


def test_genesyn_init2():
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1
    )
    eq_(
       (gs.gid, gs.taxon),
       (gs.gene_id, gs.taxon_id)
    )


def test_genesyn_calc_ratio1():
    # syn and nonsyn > 0
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=2,
        exp_syn=4,
        nonsyn=3,
        exp_nonsyn=6,
    )
    eq_(
        gs.calc_ratio(),
        1.0
    )


def test_genesyn_calc_ratio2():
    # syn = nonsyn = 0, haplotypes=True
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
    )
    eq_(
        gs.calc_ratio(haplotypes=True),
        0.0
    )


def test_genesyn_calc_ratio3a():
    # syn = nonsyn = 0
    # enough coverage
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=4
    )
    eq_(
        gs.calc_ratio(min_cov=4),
        0.0
    )


def test_genesyn_add1():
    # checks if the conserved gene_id and taxon_id is from the first instance
    gs1 = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6
    )
    gs2 = GeneSyn(
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=4
    )
    gs1.add(gs2)
    eq_(
        (gs1.gene_id, gs1.taxon_id),
        ('A', 1)
    )


def test_genesyn_add2():
    # checks if the values (syn, nonsyn, exp_syn, exp_nonsyn) are added
    gs1 = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=3,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6
    )
    gs2 = GeneSyn(
        syn=1,
        exp_syn=3,
        nonsyn=2,
        exp_nonsyn=5,
        coverage=4
    )
    gs1.add(gs2)
    eq_(
        (gs1.exp_syn, gs1.exp_nonsyn, gs1.syn, gs1.nonsyn),
        (7, 11, 4, 3)
    )


def test_genesyn_add3():
    # checks if the values if coverage is correctly added
    # coverage in gs2 is not None, in gs1 is None
    gs1 = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=3,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6
    )
    gs2 = GeneSyn(
        syn=1,
        exp_syn=3,
        nonsyn=2,
        exp_nonsyn=5,
        coverage=5
    )
    gs1.add(gs2)
    eq_(
        gs1.coverage,
        5
    )


def test_genesyn_add4():
    # checks if the values if coverage is correctly added
    # coverage in gs2 is not None, in gs1 is not None
    gs1 = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=3,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6,
        coverage=2
    )
    gs2 = GeneSyn(
        syn=1,
        exp_syn=3,
        nonsyn=2,
        exp_nonsyn=5,
        coverage=5
    )
    gs1.add(gs2)
    eq_(
        gs1.coverage,
        7
    )


def test_genesyn_add5():
    # checks if the values if coverage is correctly added
    # coverage in gs2 is None, in gs1 is None
    gs1 = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=3,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6,
    )
    gs2 = GeneSyn(
        syn=1,
        exp_syn=3,
        nonsyn=2,
        exp_nonsyn=5,
    )
    gs1.add(gs2)
    eq_(
        gs1.coverage,
        None
    )


def test_genesyn_calc_ratio3b():
    # syn = nonsyn = 0
    # coverage is None
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio3c():
    # syn = nonsyn = 0
    # coverage is < min_cov
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio(min_cov=4)))


def test_genesyn_calc_ratio4a():
    # nonsyn != 0, syn > 0
    # flag_value=True
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6,
        coverage=2
    )
    eq_(
        gs.calc_ratio(flag_value=True),
        -1
    )


def test_genesyn_calc_ratio4b():
    # nonsyn != 0, syn > 0
    # flag_value=False
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=1,
        exp_nonsyn=6,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio5a():
    # nonsyn = 0, syn != 0
    # flag_value=True
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=1,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    eq_(
        gs.calc_ratio(flag_value=True),
        -2
    )


def test_genesyn_calc_ratio5b():
    # nonsyn = 0, syn != 0
    # flag_value=False
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=1,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio6a():
    # nonsyn = syn = 0
    # flag_value=True
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    eq_(
        gs.calc_ratio(flag_value=True),
        -3
    )


def test_genesyn_calc_ratio6b():
    # nonsyn = syn = 0
    # flag_value=False
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio7a():
    # exp_nonsyn or exp_syn = 0
    # or any other case
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=4,
        nonsyn=0,
        exp_nonsyn=0,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio7b():
    # exp_nonsyn or exp_syn = 0
    # or any other case
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=0,
        nonsyn=0,
        exp_nonsyn=6,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_calc_ratio7c():
    # exp_nonsyn or exp_syn = 0
    # or any other case
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=0,
        exp_syn=0,
        nonsyn=0,
        exp_nonsyn=0,
        coverage=2
    )
    ok_(numpy.isnan(gs.calc_ratio()))


def test_genesyn_pickle_dump():
    gs = GeneSyn(
        gene_id='A',
        taxon_id=1,
        syn=1,
        exp_syn=2,
        nonsyn=3,
        exp_nonsyn=4,
        coverage=2
    )
    eq_(
        gs.__getstate__(),
        {
            'gene_id': 'A',
            'taxon_id': 1,
            'syn': 1,
            'exp_syn': 2,
            'nonsyn': 3,
            'exp_nonsyn': 4,
            'coverage': 2,
            'taxon_root': '',
        }
    )


def test_genesyn_pickle_load():
    state = {
        'gene_id': 'A',
        'taxon_id': 1,
        'syn': 1,
        'exp_syn': 2,
        'nonsyn': 3,
        'exp_nonsyn': 4,
        'coverage': 2,
        'taxon_root': '',
    }
    gs = GeneSyn()
    gs.__setstate__(state)
    eq_(
        gs.__getstate__(),
        state
    )

SNP_DATA = {
    'sample1': {
        'gene1': GeneSyn(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=6,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2
        ),  # pN/pS = 1.0
        'gene2': GeneSyn(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2
        )  # pN/pS = 2.0
    },
    'sample2': {
        'gene1': GeneSyn(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=3,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2
        ),  # pN/pS = 1.0
        'gene2': GeneSyn(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=6,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2
        ),  # pN/pS = 2.0
        'gene3': GeneSyn(
            gene_id='gene3',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2
        )  # pN/pS = 2.0
    }

}


SNP_DATA2 = {
    'sample1': {
        'gene1': GeneSyn(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=6,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2,
            coverage=4
        ),  # pN/pS = 1.0
        'gene2': GeneSyn(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=3,
            exp_nonsyn=4,
            syn=3,
            nonsyn=2,
            coverage=4
        )  # pN/pS = 2.0
    },
    'sample2': {
        'gene1': GeneSyn(
            gene_id='gene1',
            taxon_id=839,  # prevotella ruminicola
            exp_syn=3,
            exp_nonsyn=4,
            syn=4,
            nonsyn=2,
            coverage=4
        ),  # pN/pS = 1.0
        'gene2': GeneSyn(
            gene_id='gene2',
            taxon_id=838,  # prevotella genus
            exp_syn=6,
            exp_nonsyn=4,
            syn=5,
            nonsyn=2,
            coverage=3
        )  # pN/pS = 2.0
    }

}


def test_flat_sample_snps1():
    data = flat_sample_snps(SNP_DATA2, 4)
    eq_(
        (
            data['all_samples']['gene1'].syn,
            data['all_samples']['gene2'].syn
        ),
        (7, 3)
    )


def test_flat_sample_snps2():
    data = flat_sample_snps(SNP_DATA2, 3)
    eq_(
        (
            data['all_samples']['gene1'].syn,
            data['all_samples']['gene2'].syn
        ),
        (7, 8)
    )


def test_combine_sample_min_num1():
    # check if min_num is working
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    eq_(
        df.shape,
        (3, 2)
    )


def test_combine_sample_snps_min_num2():
    # check if min_num is working
    df = combine_sample_snps(SNP_DATA, 2, [lambda x: x])
    eq_(
        df.shape,
        (2, 2)
    )


def test_combine_sample_snps_values():
    # check if values are correct is working
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    eq_(
        (df.min().min(), df.max().max()),
        (0.5, 1.)
    )


def test_combine_sample_snps_index1():
    # check index_type
    df = combine_sample_snps(SNP_DATA, 1, [lambda x: x])
    eq_(
        df.index.tolist(),
        [('gene1', 839), ('gene2', 838), ('gene3', 838)]
    )


def test_combine_sample_snps_index2():
    # check index_type
    df = combine_sample_snps(
        SNP_DATA,
        1,
        [lambda x: x],
        index_type='gene'
    )
    eq_(
        df.index.tolist(),
        ['gene1', 'gene2', 'gene3']
    )


def test_combine_sample_snps_index3():
    # check index_type
    df = combine_sample_snps(
        SNP_DATA,
        1,
        [lambda x: x],
        index_type='taxon'
    )
    eq_(
        df.index.tolist(),
        [838, 839]
    )
