from nose.tools import *

from mgkit.io import gff
import misc_data


@with_setup(setup=misc_data.setup_gff_data)
def test_fromgff1():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        "KMSRIGKLPITVPAGVTVTVDENNLVTVKGPKGTLSQQVNPDITLKQEGNILTLERPTDSKPHKAMHGL",
        ann.attr['aa_seq']
    )


@with_setup(setup=misc_data.setup_gff_data)
def test_fromgff2():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(3, ann.start)


@with_setup(setup=misc_data.setup_gff_data)
def test_fromgff3():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(209, ann.end)


@with_setup(setup=misc_data.setup_gff_data)
def test_uid_fromgff_nouid1():
    # a uid is always created with fromgff, if not found, takes precedence over
    # ko_idx
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    ok_(
        len(ann.uid) != len('K02933.12503')
    )


def test_uid_fromgff_nouid2():
    # a uid is always created with fromgff, if not found, must be random, not
    # the same as another
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    ok_(
        ann.uid != '32ea1cc8-9e76-4310-8d1c-8e7890734a6b'
    )


def test_uid_fromgff_uid1():
    # a uid is not created with fromgff, if present
    line = misc_data.GFF_FILE[1]

    ann = gff.from_gff(line)

    eq_(
        ann.uid, '32ea1cc8-9e76-4310-8d1c-8e7890734a6b'
    )


def test_Annotation_dbq():
    # a dbq value is always an int
    # the same as another
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.dbq, 10
    )


def test_Annotation_get_ec1():
    # a list is returned
    line = misc_data.GFF_FILE[1]

    ann = gff.from_gff(line)

    eq_(
        ann.get_ec(), ['1.1.-', '2.2.3.4']
    )


def test_Annotation_get_ec2():
    # a level can be specified
    line = misc_data.GFF_FILE[1]

    ann = gff.from_gff(line)

    eq_(
        ann.get_ec(level=2), ['1.1', '2.2']
    )


def test_Annotation_get_ec3():
    # if no EC information is present, an empty list is returned
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.get_ec(), []
    )


def test_Annotation_get_mapping1():
    # a list is returned
    line = misc_data.GFF_FILE[2]

    ann = gff.from_gff(line)

    eq_(
        ann.get_mapping('test'), ['12345']
    )


def test_Annotation_get_mapping2():
    # an empty list is returned if not mapping is found
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.get_mapping('test'), []
    )


def test_Annotation_add_exp_syn_count():
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    ann.add_exp_syn_count(misc_data.NUC_SEQS['contig-1327918'])

    eq_(
        (141, 480),
        (ann.exp_syn, ann.exp_nonsyn)
    )


def test_Annotation_add_gc_content():
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    ann.add_gc_content(misc_data.NUC_SEQS['contig-1327918'])

    eq_(
        0.5314009661835749,
        ann.get_attr('gc_cont', float)
    )


def test_Annotation_add_gc_ratio():
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    ann.add_gc_ratio(misc_data.NUC_SEQS['contig-1327918'])

    eq_(
        0.8818181818181818,
        ann.get_attr('gc_ratio', float)
    )


def test_Annotation_to_gff():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.attr, gff.from_gff(ann.to_gff()).attr
    )


def test_Annotation_to_gtf():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.gene_id, gff.from_gff(ann.to_gtf()).attr['transcript_id']
    )


def test_Annotation_sample_coverage():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        int(ann.attr['t1_b3_cov']), ann.sample_coverage['t1_b3']
    )


def test_Annotation_get_number_of_samples1():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.get_number_of_samples(min_cov=0), 14
    )


def test_Annotation_get_number_of_samples2():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.get_number_of_samples(min_cov=15), 8
    )


@with_setup(setup=misc_data.setup_aaseq_data)
@with_setup(setup=misc_data.setup_hmmer_data)
def test_from_hmmer():
    checks = (
        'contig-1442648',
        'K00001_4479_poaceae',
        693,
        894,
        'K00001',
        4479,
        'poaceae'
    )
    ann = gff.from_hmmer(
        misc_data.HMMER_FILE[0], misc_data.AA_SEQS
    )

    eq_(
        (
            ann.seq_id,
            ann.attr['name'],
            ann.attr['aa_from'],
            ann.attr['aa_to'],
            ann.gene_id,
            ann.taxon_id,
            ann.attr['taxon_name']
        ),
        checks
    )


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
