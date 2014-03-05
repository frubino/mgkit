from nose.tools import *

from mgkit.snps.classes import GeneSyn

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
    #syn and nonsyn > 0
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
    #syn = nonsyn = 0, haplotypes=True
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
    #syn = nonsyn = 0
    #enough coverage
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


def test_genesyn_calc_ratio3b():
    #syn = nonsyn = 0
    #coverage is None
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
    #syn = nonsyn = 0
    #coverage is < min_cov
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
    #nonsyn != 0, syn > 0
    #flag_value=True
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
    #nonsyn != 0, syn > 0
    #flas_value=False
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
    #nonsyn = 0, syn != 0
    #flag_value=True
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
    #nonsyn = 0, syn != 0
    #flas_value=False
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
    #nonsyn = syn = 0
    #flag_value=True
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
    #nonsyn = syn = 0
    #flas_value=False
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
    #exp_nonsyn or exp_syn = 0
    #or any other case
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
    #exp_nonsyn or exp_syn = 0
    #or any other case
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
    #exp_nonsyn or exp_syn = 0
    #or any other case
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
