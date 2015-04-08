from nose.tools import *

from mgkit.io import gff
from mgkit.utils import sequence
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
        ann.get_ec(), set(['1.1.-', '2.2.3.4'])
    )


def test_Annotation_get_ec2():
    # a level can be specified
    line = misc_data.GFF_FILE[1]

    ann = gff.from_gff(line)

    eq_(
        ann.get_ec(level=2), set(['1.1', '2.2'])
    )


def test_Annotation_get_ec2__duplicates():
    # a level can be specified
    line = misc_data.GFF_FILE[1]

    ann = gff.from_gff(line)
    ann.attr['EC'] = ann.attr['EC'] + ',2.2.3.1'

    eq_(
        ann.get_ec(level=2), set(['1.1', '2.2'])
    )


def test_Annotation_get_ec3():
    # if no EC information is present, an empty list is returned
    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.get_ec(), set([])
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


def test_Annotation_to_gtf1():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.uid, gff.from_gff(ann.to_gtf()).attr['transcript_id']
    )


def test_Annotation_to_gtf2():

    line = misc_data.GFF_FILE[0]

    ann = gff.from_gff(line)

    eq_(
        ann.gene_id, gff.from_gff(ann.to_gtf(gene_id_attr='ko')).attr['transcript_id']
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


def test_Annotation_get_nuc_seq1():
    ann = gff.Annotation(start=1, end=40, strand='+')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq),
        seq
    )


def test_Annotation_get_nuc_seq2():
    ann = gff.Annotation(start=2, end=40, strand='+')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq),
        seq[1:]
    )


def test_Annotation_get_nuc_seq3():
    ann = gff.Annotation(start=2, end=39, strand='+')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq),
        seq[1:-1]
    )


def test_Annotation_get_nuc_seq_reverse1():
    ann = gff.Annotation(start=1, end=40, strand='-')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq),
        seq
    )


def test_Annotation_get_nuc_seq_reverse2():
    ann = gff.Annotation(start=2, end=39, strand='-')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq, reverse=True),
        sequence.reverse_complement(seq[1:-1])
    )


def test_Annotation_get_nuc_seq_snp():
    ann = gff.Annotation(start=2, end=39, strand='+')
    seq = 'ACTG' * 10

    eq_(
        ann.get_nuc_seq(seq, snp=(2, 'C')),
        sequence.get_variant_sequence(seq[1:-1], (2, 'C'))
    )


def test_Annotation_get_aa_seq():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.get_aa_seq(seq),
        sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=0,
            reverse=False
        )
    )


def test_Annotation_get_aa_seq__start1():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.get_aa_seq(seq, start=None),
        sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=ann.phase,
            reverse=False
        )
    )


def test_Annotation_get_aa_seq__start2():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    eq_(
        ann.get_aa_seq(seq, start=None),
        sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=ann.phase,
            reverse=False
        )
    )


def test_Annotation_get_aa_seq__snp():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    eq_(
        ann.get_aa_seq(seq, start=None, snp=(1, 'C')),
        sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True, snp=(1, 'C')),
            start=ann.phase,
            reverse=False
        )
    )


def test_Annotation_is_syn1_3():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 3, 'C'),
        True
    )


def test_Annotation_is_syn1_1():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 1, 'C'),
        False
    )


def test_Annotation_is_syn2_1():
    #second position on reference, first base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 2, 'A'),
        False
    )


def test_Annotation_is_syn2_2():
    #second position on reference, second base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 3, 'A'),
        False
    )


def test_Annotation_is_syn2_3():
    #second position on reference, third base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 4, 'A'),
        True
    )


def test_Annotation_is_syn3_1():
    #first position on reference, second codon, first base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 4, 'A'),
        False
    )


def test_Annotation_is_syn3_2():
    #first position on reference, second codon, second base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 5, 'T'),
        False
    )


def test_Annotation_is_syn3_3():
    #first position on reference, second codon, third base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 6, 'T'),
        True
    )


def test_Annotation_is_syn4_1():
    #first position on reference, third codon, first base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 7, 'A'),
        False
    )


def test_Annotation_is_syn4_2():
    #first position on reference, third codon, second base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 8, 'T'),
        False
    )


def test_Annotation_is_syn4_3():
    #first position on reference, third codon, third base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 9, 'T'),
        False
    )


def test_Annotation_is_syn__start1_1():
    #first position on reference, second codon, first base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 5, 'G', start=1),
        False
    )


def test_Annotation_is_syn__start1_2():
    #first position on reference, second codon, second base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 6, 'A', start=1),
        False
    )


def test_Annotation_is_syn__start1_3():
    #first position on reference, second codon, third base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 6, 'A', start=1),
        False
    )


def test_Annotation_is_syn__start2_1():
    #first position on reference, second codon, first base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 6, 'A', start=None),
        False
    )


def test_Annotation_is_syn__start3_1():
    #first position on reference, third codon, third base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 10, 'T', start=None),
        True
    )


def test_Annotation_is_syn__start3_2():
    #first position on reference, third codon, second base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    eq_(
        ann.is_syn(seq, 9, 'C', start=None),
        False
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


def test_genomicrange_contains1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        1 in gen_range1,
        False
    )


def test_genomicrange_contains2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        20 in gen_range1,
        True
    )


def test_genomicrange_contains3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        30 in gen_range1,
        True
    )


def test_genomicrange_contains4():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        31 in gen_range1,
        False
    )


def test_genomicrange_get_relative_pos1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(20),
        1
    )


def test_genomicrange_get_relative_pos2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(30),
        11
    )


def test_genomicrange_get_relative_pos3():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    eq_(
        gen_range1.get_relative_pos(25),
        6
    )


@raises(ValueError)
def test_genomicrange_get_relative_pos_fail1():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    gen_range1.get_relative_pos(19)


@raises(ValueError)
def test_genomicrange_get_relative_pos_fail2():
    gen_range1 = gff.GenomicRange(start=20, end=30)

    gen_range1.get_relative_pos(31)


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
