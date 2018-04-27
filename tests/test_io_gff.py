import pytest

from mgkit.io import gff, fasta

@pytest.fixture
def nucseq(shared_datadir):
    return dict(
        fasta.load_fasta((shared_datadir / 'test-seq-nuc.fa').open())
    )

@pytest.fixture
def gff_file(shared_datadir):
    return (shared_datadir / 'test.gff').open().readlines()


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

    assert (147, 474) == (ann.exp_syn, ann.exp_nonsyn)


def test_Annotation_add_gc_content(gff_file, nucseq):

    ann = gff.from_gff(gff_file[0])
    ann.add_gc_content(nucseq['contig-1327918'])

    assert 0.5314009661835749 == ann.get_attr('gc_cont', float)


def test_Annotation_add_gc_ratio(gff_file, nucseq):

    ann = gff.from_gff(gff_file[0])
    ann.add_gc_ratio(nucseq['contig-1327918'])

    assert 0.8818181818181818 == ann.get_attr('gc_ratio', float)
