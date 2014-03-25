from nose.tools import *

import misc_data

from mgkit.io import gff


def test_gffattributesdict_init():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann2 = gff.GFFAttributesDict()
    ann2.ko_idx = 'test'
    ann2.cov = 3
    eq_(ann1, ann2)


def test_gffattributesdict_setattr():
    ann1 = gff.GFFAttributesDict()
    ann2 = gff.GFFAttributesDict()
    ann1['ko_idx'] = 'test'
    ann1['cov'] = 3
    ann2.ko_idx = 'test'
    ann2.cov = 3
    eq_(ann1, ann2)


def test_gffattributesdict_getattr():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann2 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    eq_(ann1['ko_idx'], ann2.ko_idx)


def test_gffattributesdict_hash():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann2 = gff.GFFAttributesDict()
    ann2.ko_idx = 'test'
    ann2.cov = 3
    eq_(hash(ann1), hash(ann2))


def test_gffattributesdict_hash2():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann2 = gff.GFFAttributesDict()
    ann2.ko_idx = 'test'
    ann2.cov = 3
    ann1.calc_hash()
    ann2.calc_hash()
    eq_(ann1._hash, ann2._hash)


def test_gffattributesdict_hash3():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann2 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    ann1.calc_hash()
    ann2.calc_hash()
    ann2['cov'] = 9
    eq_(ann1._hash, ann2._hash)


def test_gffattributesdict_to_string():
    ann1 = gff.GFFAttributesDict(ko_idx='test', cov=3)
    eq_(ann1.to_string(), 'cov="3";ko_idx="test"')


@with_setup(setup=misc_data.setup_gff_data)
def test_basegffdict_parse_line():

    line = misc_data.GFF_FILE[0]

    ann = gff.BaseGFFDict(line)

    eq_(
        "KMSRIGKLPITVPAGVTVTVDENNLVTVKGPKGTLSQQVNPDITLKQEGNILTLERPTDSKPHKAMHGL",
        ann.attributes.aa_seq
    )


@with_setup(setup=misc_data.setup_gff_data)
def test_basegffdict_parse_line2():

    line = misc_data.GFF_FILE[0]

    ann = gff.BaseGFFDict(line)

    eq_(209, ann.feat_to)


@with_setup(setup=misc_data.setup_gff_data)
def test_basegffdict_calc_hash():

    line = misc_data.GFF_FILE[0]

    ann1 = gff.BaseGFFDict(line)
    ann2 = gff.BaseGFFDict(line)

    eq_(hash(ann1), hash(ann2))


@with_setup(setup=misc_data.setup_gff_data)
def test_basegffdict_calc_hash2():

    line1 = misc_data.GFF_FILE[0]
    line2 = misc_data.GFF_FILE[1]

    ann1 = gff.BaseGFFDict(line1)
    ann2 = gff.BaseGFFDict(line2)

    assert hash(ann1) != hash(ann2)


@with_setup(setup=misc_data.setup_gff_data)
def test_basegffdict_to_string():

    line = misc_data.GFF_FILE[0]

    ann1 = gff.BaseGFFDict(line)
    ann2 = gff.BaseGFFDict(ann1.to_string())

    eq_(hash(ann1), hash(ann2))


@with_setup(setup=misc_data.setup_nucseq_data)
@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_gffkegg_from_hmmer():
    checks = (
        'contig-1442648',
        'K00001_4479_poaceae',
        693,
        894,
        'K00001',
        '4479',
        'poaceae'
    )
    ann = gff.GFFKegg.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS, misc_data.NUC_SEQS
    )

    eq_(
        (
            ann.seq_id,
            ann.attributes.name,
            ann.attributes.aa_from,
            ann.attributes.aa_to,
            ann.attributes.ko,
            ann.attributes.taxon_id,
            ann.attributes.taxon
        ),
        checks
    )


@with_setup(setup=misc_data.setup_nucseq_data)
@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_gffkegg_to_gff():
    ann = gff.GFFKegg.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS, misc_data.NUC_SEQS
    )
    ann.attributes.ko_idx = 'K00001.1'
    ann = ann.to_gtf()

    eq_(
        (ann.attributes.gene_id, ann.attributes.transcript_id),
        ('K00001.1', 'K00001.1')
    )


@with_setup(setup=misc_data.setup_nucseq_data)
@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_gffkegg_get_taxon_id1():
    ann = gff.GFFKegg.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS, misc_data.NUC_SEQS
    )
    ann.attributes.taxon_id = 12
    ann.attributes.blast_taxon_idx = 1
    eq_(ann.get_taxon_id(), 1)


@with_setup(setup=misc_data.setup_nucseq_data)
@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_gffkegg_get_taxon_id2():
    ann = gff.GFFKegg.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS, misc_data.NUC_SEQS
    )
    ann.attributes.taxon_id = 12
    ann.attributes.blast_taxon_idx = 1
    eq_(ann.get_taxon_id(prefer_blast=False), 12)


def test_gff_glimmer3_line1():
    header = 'sequence0001'
    line = 'orf00001       66      611  +3     6.08'
    annotation = gff.BaseGFFDict.from_glimmer3(header, line)

    eq_(
        (
            annotation.seq_id,
            annotation.feat_from,
            annotation.feat_to,
            annotation.score,
            annotation.strand,
            annotation.phase,
            annotation.attributes.ID,
            annotation.attributes.frame
        ),
        (
            header,
            66,
            611,
            6.08,
            '+',
            2,
            'orf00001',
            '+3'
        )
    )


def test_gff_glimmer3_line2():
    header = 'sequence0001'
    line = 'orf00001       66      11  -2     6.08'
    annotation = gff.BaseGFFDict.from_glimmer3(header, line)

    eq_(
        (
            annotation.feat_from,
            annotation.feat_to,
            annotation.strand,
            annotation.phase,
            annotation.attributes.frame
        ),
        (
            11,
            66,
            '-',
            1,
            '-2'
        )
    )


@with_setup(setup=misc_data.setup_glimmer3_data)
def test_gff_glimmer_parse():
    annotations = list(gff.parse_glimmer3(misc_data.GLIMMER3_FILE))

    eq_(
        (
            sum(1 for annotation in annotations if annotation.seq_id == 'scaff01'),
            sum(1 for annotation in annotations if annotation.seq_id == 'scaff21'),
            annotations[4].feat_from,
            annotations[6].feat_from,
            annotations[7].feat_from,
            annotations[8].feat_from,
        ),
        (
            7,
            2,
            3340,
            6079,
            76,
            60
        )
    )
