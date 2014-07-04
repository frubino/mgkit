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
    annotation = gff.from_glimmer3(header, line)

    eq_(
        (
            annotation.seq_id,
            annotation.start,
            annotation.end,
            annotation.score,
            annotation.strand,
            annotation.phase,
            annotation.attr['orf_id'],
            annotation.attr['frame']
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
    annotation = gff.from_glimmer3(header, line)

    eq_(
        (
            annotation.start,
            annotation.end,
            annotation.strand,
            annotation.phase,
            annotation.attr['frame']
        ),
        (
            11,
            66,
            '-',
            1,
            '-2'
        )
    )


@with_setup(setup=misc_data.setup_nucseq_data)
def test_gff_from_sequence1():
    annotation = gff.from_sequence(
        'contig-110637',
        misc_data.NUC_SEQS['contig-110637']
    )
    eq_(
        (
            annotation.seq_id,
            annotation.start,
            annotation.end,
        ),
        (
            'contig-110637',
            1,
            len(misc_data.NUC_SEQS['contig-110637'])
        )
    )


def test_genomicrange_union1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 20)
    )


def test_genomicrange_union3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 20
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union4():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 20
    gen_range2.end = 30

    gen_range_u = gen_range2.union(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 30)
    )


def test_genomicrange_union_fail1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq2'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_union_fail2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '-'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_union_fail3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 21
    gen_range2.end = 30

    gen_range_u = gen_range1.union(gen_range2)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 19
    gen_range2.end = 30

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (19, 20)
    )


def test_genomicrange_intersect2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 15
    gen_range2.end = 30

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (15, 20)
    )


def test_genomicrange_intersect3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 15
    gen_range2.end = 30

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (15, 20)
    )


def test_genomicrange_intersect4():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        (gen_range_u.start, gen_range_u.end),
        (10, 20)
    )


def test_genomicrange_intersect5():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        len(gen_range_u),
        len(gen_range2)
    )


def test_genomicrange_intersect6():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range1.intersect(gen_range2)

    eq_(
        len(gen_range_u),
        len(gen_range2)
    )


def test_genomicrange_intersect_fail1():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq2'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect_fail2():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '-'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 10
    gen_range2.end = 20

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def test_genomicrange_intersect_fail3():
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = 30
    gen_range2.end = 40

    gen_range_u = gen_range2.intersect(gen_range1)

    eq_(
        gen_range_u,
        None
    )


def elongate_data():
    seq_id = 'test1'

    test_ann = [
        gff.GenomicRange(seq_id=seq_id, start=1, end=10),
        gff.GenomicRange(seq_id=seq_id, start=10, end=15),
        gff.GenomicRange(seq_id=seq_id, start=16, end=18),
    ]
    return test_ann


def test_annotation1():
    ann = gff.Annotation(start=1, end=10)
    eq_(len(ann), 10)


def test_genomicrange_misc1():
    ann = gff.GenomicRange(start=1, end=10)
    eq_(len(ann), 10)
