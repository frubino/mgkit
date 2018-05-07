import pytest

from mgkit.io import gff, fasta
from mgkit.io import open_file
from mgkit.utils import sequence

@pytest.fixture
def nucseq(shared_datadir):
    return dict(
        fasta.load_fasta(str(shared_datadir / 'test-seq-nuc.fa'))
    )


@pytest.fixture
def aaseq(shared_datadir):
    return dict(
        fasta.load_fasta(str(shared_datadir / 'test-seq-aa.fa'))
    )


@pytest.fixture
def gff_file(shared_datadir):
    return open_file(str(shared_datadir / 'test.gff'), 'rb').readlines()


@pytest.fixture
def hmmer_file(shared_datadir):
    return open_file(str(shared_datadir / 'test-hmmer-dom.txt'), 'rb').readlines()


@pytest.fixture
def glimmer_file(shared_datadir):
    return open_file(str(shared_datadir / 'glimmer3.txt'), 'rb').readlines()


@pytest.fixture
def keggmod_file(shared_datadir):
    return open_file(str(shared_datadir / 'kmod-entry1.txt'), 'rb').readlines()


def test_fromgff1(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert "KMSRIGKLPITVPAGVTVTVDENNLVTVKGPKGTLSQQVNPDITLKQEGNILTLERPTDSKPHKAMHGL" ==  ann.attr['aa_seq']


def test_fromgff2(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert 3 == ann.start


def test_fromgff3(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert 209 == ann.end


def test_fromgff_nouid1(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert len(ann.uid) != len('K02933.12503')


def test_fromgff_nouid2(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.uid != '32ea1cc8-9e76-4310-8d1c-8e7890734a6b'


def test_fromgff_uid1(gff_file):

    ann = gff.from_gff(gff_file[1])

    assert ann.uid == '32ea1cc8-9e76-4310-8d1c-8e7890734a6b'


def test_Annotation_dbq(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.dbq == 10


def test_Annotation_ec1(gff_file):

    ann = gff.from_gff(gff_file[1])

    assert ann.get_ec() == set(['1.1.-', '2.2.3.4'])


def test_Annotation_ec2(gff_file):

    ann = gff.from_gff(gff_file[1])

    assert ann.get_ec(level=2) == set(['1.1', '2.2'])


def test_Annotation_ec2_dup(gff_file):

    ann = gff.from_gff(gff_file[1])
    ann.attr['EC'] = ann.attr['EC'] + ',2.2.3.1'

    assert ann.get_ec(level=2) == set(['1.1', '2.2'])


def test_Annotation_ec3(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.get_ec() == set()


def test_Annotation_get_mapping1(gff_file):

    ann = gff.from_gff(gff_file[2])

    assert ann.get_mapping('test') == ['12345']


def test_Annotation_get_mapping2(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.get_mapping('test') == []


def test_Annotation_add_exp_syn_count(gff_file, nucseq):

    ann = gff.from_gff(gff_file[0])
    ann.add_exp_syn_count(nucseq['contig-1327918'])

    assert (141, 480) == (ann.exp_syn, ann.exp_nonsyn)


def test_Annotation_add_gc_content(gff_file, nucseq):

    ann = gff.from_gff(gff_file[0])
    ann.add_gc_content(nucseq['contig-1327918'])

    assert 0.5314009661835749 == ann.get_attr('gc_cont', float)


def test_Annotation_add_gc_ratio(gff_file, nucseq):

    ann = gff.from_gff(gff_file[0])
    ann.add_gc_ratio(nucseq['contig-1327918'])

    assert 0.8818181818181818 == ann.get_attr('gc_ratio', float)


def test_Annotation_to_gff(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert gff.from_gff(ann.to_gff()).attr == ann.attr


def test_Annotation_to_gtf1(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert gff.from_gff(ann.to_gtf()).attr['transcript_id'] == ann.uid


def test_Annotation_to_gtf2(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert gff.from_gff(ann.to_gtf(gene_id_attr='ko')).attr['transcript_id'] == ann.gene_id


def test_Annotation_sample_coverage(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert int(ann.attr['t1_b3_cov']) == ann.sample_coverage['t1_b3']


def test_Annotation_get_number_of_samples1(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.get_number_of_samples(min_cov=0) == 14


def test_Annotation_get_number_of_samples2(gff_file):

    ann = gff.from_gff(gff_file[0])

    assert ann.get_number_of_samples(min_cov=15) == 8


def test_Annotation_get_nuc_seq1():
    ann = gff.Annotation(start=1, end=40, strand='+')
    seq = 'ACTG' * 10

    assert ann.get_nuc_seq(seq) == seq


def test_Annotation_get_nuc_seq2():
    ann = gff.Annotation(start=2, end=40, strand='+')
    seq = 'ACTG' * 10

    assert ann.get_nuc_seq(seq) == seq[1:]


def test_Annotation_get_nuc_seq3():
    ann = gff.Annotation(start=2, end=39, strand='+')
    seq = 'ACTG' * 10

    assert ann.get_nuc_seq(seq) == seq[1:-1]


def test_Annotation_get_nuc_seq_reverse1():
    ann = gff.Annotation(start=1, end=40, strand='-')
    seq = 'ACTG' * 10

    assert ann.get_nuc_seq(seq) == seq


def test_Annotation_get_nuc_seq_reverse2():
    ann = gff.Annotation(start=2, end=39, strand='-')
    seq = 'ACTG' * 10

    assert ann.get_nuc_seq(seq, reverse=True) == sequence.reverse_complement(seq[1:-1])


def test_Annotation_get_nuc_seq_snp():
    ann = gff.Annotation(start=2, end=39, strand='+')
    seq = 'ACTG' * 10

    ann.get_nuc_seq(seq, snp=(2, 'C')) == sequence.get_variant_sequence(seq[1:-1], (2, 'C'))


def test_Annotation_get_aa_seq():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert ann.get_aa_seq(seq) == sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=0,
            reverse=False
        )


def test_Annotation_get_aa_seq__start1():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert ann.get_aa_seq(seq, start=None) == sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=ann.phase,
            reverse=False
        )


def test_Annotation_get_aa_seq__start2():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    assert ann.get_aa_seq(seq, start=None) == sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True),
            start=ann.phase,
            reverse=False
        )


def test_Annotation_get_aa_seq__snp():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    assert ann.get_aa_seq(seq, start=None, snp=(1, 'C')) == sequence.translate_sequence(
            ann.get_nuc_seq(seq, reverse=True, snp=(1, 'C')),
            start=ann.phase,
            reverse=False
        )


def test_Annotation_is_syn1_3():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert ann.is_syn(seq, 3, 'C')


def test_Annotation_is_syn1_1():
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 1, 'C')


def test_Annotation_is_syn2_1():
    # second position on reference, first base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 2, 'A')


def test_Annotation_is_syn2_2():
    # second position on reference, second base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 3, 'A')


def test_Annotation_is_syn2_3():
    # second position on reference, third base in codon
    ann = gff.Annotation(start=2, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert ann.is_syn(seq, 4, 'A')


def test_Annotation_is_syn3_1():
    # first position on reference, second codon, first base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 4, 'A')


def test_Annotation_is_syn3_2():
    # first position on reference, second codon, second base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 5, 'T')


def test_Annotation_is_syn3_3():
    # first position on reference, second codon, third base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert ann.is_syn(seq, 6, 'T')


def test_Annotation_is_syn4_1():
    # first position on reference, third codon, first base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 7, 'A')


def test_Annotation_is_syn4_2():
    # first position on reference, third codon, second base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 8, 'T')


def test_Annotation_is_syn4_3():
    # first position on reference, third codon, third base in codon
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 9, 'T')


def test_Annotation_is_syn__start1_1():
    # first position on reference, second codon, first base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 5, 'G', start=1)


def test_Annotation_is_syn__start1_2():
    # first position on reference, second codon, second base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 6, 'A', start=1)


def test_Annotation_is_syn__start1_3():
    # first position on reference, second codon, third base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=0)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 6, 'A', start=1)


def test_Annotation_is_syn__start2_1():
    # first position on reference, second codon, first base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 6, 'A', start=None)


def test_Annotation_is_syn__start3_1():
    # first position on reference, third codon, third base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    assert ann.is_syn(seq, 10, 'T', start=None)


def test_Annotation_is_syn__start3_2():
    # first position on reference, third codon, second base in codon, phase=1
    ann = gff.Annotation(start=1, end=40, strand='+', phase=1)
    seq = 'ACTG' * 10

    assert not ann.is_syn(seq, 9, 'C', start=None)


def test_from_hmmer(hmmer_file, aaseq):
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
        hmmer_file[0], aaseq
    )

    assert (
            ann.seq_id,
            ann.attr['name'],
            ann.attr['aa_from'],
            ann.attr['aa_to'],
            ann.gene_id,
            ann.taxon_id,
            ann.attr['taxon_name']
        ) == checks


def test_gff_glimmer3_line1():
    header = 'sequence0001'
    line = 'orf00001       66      611  +3     6.08'
    annotation = gff.from_glimmer3(header, line)

    assert (
            annotation.seq_id,
            annotation.start,
            annotation.end,
            annotation.score,
            annotation.strand,
            annotation.phase,
            annotation.attr['orf_id'],
            annotation.attr['frame']
        ) == (
            header,
            66,
            611,
            6.08,
            '+',
            2,
            'orf00001',
            '+3'
        )


def test_gff_glimmer3_line2():
    header = 'sequence0001'
    line = 'orf00001       66      11  -2     6.08'
    annotation = gff.from_glimmer3(header, line)

    assert (
            annotation.start,
            annotation.end,
            annotation.strand,
            annotation.phase,
            annotation.attr['frame']
        ) == (
            11,
            66,
            '-',
            1,
            '-2'
        )


def test_gff_from_sequence1(nucseq):
    annotation = gff.from_sequence(
        'contig-110637',
        nucseq['contig-110637']
    )
    assert (
            annotation.seq_id,
            annotation.start,
            annotation.end,
        ) == (
            'contig-110637',
            1,
            len(nucseq['contig-110637'])
        )


@pytest.mark.parametrize(
    "start,end,coords,result",
    [
        (20, 30, 1, False),
        (20, 30, 20, True),
        (20, 30, 30, True),
        (20, 30, 25, True),
        (20, 30, (25, 30), True),
        (20, 30, (25, 31), False),
        (20, 30, (19, 30), False),
        (20, 30, (30, 25), True),
    ]
)
def test_genomicrange_contains(start, end, coords, result):
    gen_range1 = gff.GenomicRange(start=start, end=end)

    assert (coords in gen_range1) == result


@pytest.mark.parametrize(
    "coords1,coords2,result",
    [
        ((20, 30), (25, 30), True),
        ((20, 30), (30, 25), True),
        ((20, 30), (19, 30), False),
        ((20, 30), (25, 31), False),
    ]
)
def test_genomicrange_contains(coords1, coords2, result):
    gen_range1 = gff.GenomicRange(start=coords1[0], end=coords1[1])
    gen_range2 = gff.GenomicRange(start=coords2[0], end=coords2[1])

    assert (gen_range2 in gen_range1) == result


@pytest.mark.parametrize(
    "coords1,coords2,result",
    [
        ((20, 30), (25, 30), True),
        ((20, 30), (30, 25), True),
        ((20, 30), (19, 30), False),
        ((20, 30), (25, 31), False),
    ]
)
def test_genomicrange_contains_annotation1(coords1, coords2, result):
    gen_range1 = gff.GenomicRange(start=coords1[0], end=coords1[1])

    assert (gff.Annotation(start=coords2[0], end=coords2[1]) in gen_range1) == result
