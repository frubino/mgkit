import pytest
import pathlib
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


def test_Annotation_to_file(gff_file, tmpdir):

    ann = gff.from_gff(gff_file[0])

    file_name = (tmpdir / 'test-write.gff').strpath
    file_handle = open_file(file_name, 'wb')
    ann.to_file(file_handle)
    file_handle.close()

    ann2 = next(gff.parse_gff(file_name))

    assert ann == ann2


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


def test_split_gff_file1(tmpdir, shared_datadir):
    gff_file = str(shared_datadir / 'test.gff')

    gff.split_gff_file(gff_file, (tmpdir / 'test{}.gff').strpath, 2)

    files = list(pathlib.Path(tmpdir.strpath).glob('*.gff'))

    assert len(files) == 2


def test_split_gff_file2(tmpdir, shared_datadir):

    gff_file = str(shared_datadir / 'test.gff')

    gff.split_gff_file(gff_file, (tmpdir / 'test{}.gff').strpath, 2)

    files = list(
        str(path)
        for path in pathlib.Path(tmpdir.strpath).glob('*.gff')
    )

    count1 = sum(1 for x in gff.parse_gff(gff_file))
    count2 = sum(1 for x in gff.parse_gff_files(files))

    assert count1 == count2


@pytest.mark.parametrize(
    "start1,end1,start2,end2,result",
    [
        (20, 30, 25, 30, True),
        (20, 30, 30, 25, True),
        (20, 30, 19, 30, False),
        (20, 30, 25, 31, False),
    ]
)
def test_genomicrange_contains_annotation(start1, end1, start2, end2, result):
    gen_range1 = gff.GenomicRange(start=start1, end=end1)

    assert (gff.Annotation(start=start2, end=end2) in gen_range1) == result


@pytest.mark.parametrize(
    "start,end,pos,relpos",
    [
        (20, 30, 20, 1),
        (20, 30, 30, 11),
        (20, 30, 25, 6),
    ]
)
def test_genomicrange_get_relative_pos(start, end, pos, relpos):
    gen_range1 = gff.GenomicRange(start=start, end=end)

    gen_range1.get_relative_pos(pos) == relpos


@pytest.mark.parametrize(
    "start,end,pos",
    [
        (20, 30, 19),
        (20, 30, 31),
    ]
)
def test_genomicrange_get_relative_pos_fail(start, end, pos):
    gen_range1 = gff.GenomicRange(start=start, end=end)

    with pytest.raises(ValueError):
        gen_range1.get_relative_pos(pos)


@pytest.mark.parametrize(
    "start2,end2,ustart,uend",
    [
        (19, 30, 10, 30),
        (10, 20, 10, 20),
        (20, 30, 10, 30),
        (21, 30, 10, 30),
    ]
)
def test_genomicrange_union(start2, end2, ustart, uend):
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = 'seq1'
    gen_range2.strand = '+'
    gen_range2.start = start2
    gen_range2.end = end2

    gen_range_u = gen_range1.union(gen_range2)

    assert (gen_range_u.start, gen_range_u.end) == (ustart, uend)


@pytest.mark.parametrize(
    "start2,end2,seq2,strand2",
    [
        (19, 30, 'seq2', '+'),
        (19, 30, 'seq1', '-'),
        (21, 30, 'seq1', '-'),
    ]
)
def test_genomicrange_union_fail(start2, end2, seq2, strand2):
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = seq2
    gen_range2.strand = strand2
    gen_range2.start = start2
    gen_range2.end = end2

    gen_range_u = gen_range1.union(gen_range2)

    assert gen_range_u is None


@pytest.mark.parametrize(
    "start2,end2,seq2,strand2,result",
    [
        (19, 30, 'seq1', '+', (19, 20)),
        (15, 30, 'seq1', '+', (15, 20)),
        (10, 20, 'seq1', '+', (10, 20)),
    ]
)
def test_genomicrange_intersect1(start2, end2, seq2, strand2, result):
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = 'seq1'
    gen_range1.strand = '+'
    gen_range1.start = 10
    gen_range1.end = 20
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = seq2
    gen_range2.strand = strand2
    gen_range2.start = start2
    gen_range2.end = end2

    gen_range_u = gen_range1.intersect(gen_range2)

    assert (gen_range_u.start, gen_range_u.end) == result



def test_genomicrange_intersect2():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range2.intersect(gen_range1)

    assert len(gen_range_u) == len(gen_range2)


def test_genomicrange_intersect3():
    gen_range1 = gff.GenomicRange(seq_id='seq1', strand='+', start=10, end=20)
    gen_range2 = gff.GenomicRange(seq_id='seq1', strand='+', start=12, end=18)

    gen_range_u = gen_range1.intersect(gen_range2)

    assert len(gen_range_u) == len(gen_range2)


@pytest.mark.parametrize(
    "seq1,strand1,s1,e1,seq2,strand2,s2,e2",
    [
        ('seq2', '+', 10, 20, 'seq1', '+', 10, 20),
        ('seq1', '-', 10, 20, 'seq1', '+', 10, 20),
        ('seq1', '+', 10, 20, 'seq1', '+', 30, 40),
    ]
)
def test_genomicrange_intersect_fail(seq1, strand1, s1, e1, seq2, strand2, s2, e2):
    gen_range1 = gff.GenomicRange()
    gen_range1.seq_id = seq1
    gen_range1.strand = strand1
    gen_range1.start = s1
    gen_range1.end = e1
    gen_range2 = gff.GenomicRange()
    gen_range2.seq_id = seq2
    gen_range2.strand = strand2
    gen_range2.start = s2
    gen_range2.end = e2

    gen_range_u = gen_range2.intersect(gen_range1)

    assert gen_range_u is None


@pytest.fixture
def elongate_data():
    seq_id = 'test1'

    test_ann = [
        gff.Annotation(seq_id=seq_id, start=1, end=10),
        gff.Annotation(seq_id=seq_id, start=10, end=15),
        gff.Annotation(seq_id=seq_id, start=16, end=18),
    ]
    return test_ann


def test_elongate_annotations(elongate_data):
    assert gff.elongate_annotations(elongate_data) == {(1, 18)}


def test_annotation_length1():
    ann = gff.Annotation(start=1, end=10)
    assert len(ann) == 10


def test_annotation_length2():
    ann = gff.Annotation(start=1, end=10)
    assert ann.length == 10


def test_genomicrange_length1():
    ann = gff.GenomicRange(start=1, end=10)
    assert len(ann) == 10


def test_genomicrange_length2():
    ann = gff.GenomicRange(start=1, end=10)
    with pytest.raises(AttributeError):
        ann.length
