import pytest
from mgkit.io import fasta
from mgkit.io.utils import open_file
import pathlib

@pytest.fixture
def nucseq(shared_datadir):
    return str(shared_datadir / 'test-seq-nuc.fa')


@pytest.fixture
def oneseq(shared_datadir):
    return str(shared_datadir / 'one-seq.fa')


def test_load_fasta1(nucseq):
    seq_id, seq = next(fasta.load_fasta(nucseq))
    assert seq == u'ATCGATGTCTGCCGCAATGACGGTGGCACCACGCTGGGCAGAGCGTACGACGGCACCCATACCGACCATGCCACAGCCAATCACCATGACGACATCGATATCTGTTACCTGGGCACGACTGACAGCATGGAAACCCACACTCATCGGTTCGATTAATGCACAGGTGCGGGGTGTCAGCAGTCCTGCGGGGATGACTTTCTCCCAGGGCAGGGCGAGATACTCACACATGGCTCCCCAGCGCTGCACACCTAATGTCTGGTTGTGCTCGCAGGCATTGACACGGTCGTTACGGCATGACGCACATTTTCCACAGTTGGTGTAGGGGTTGACGGTGACGGTCATACCGGGCTTCAGTCCCTCAGGTACGTTCTTGCCAATCTTGACAATCTCCGCACCTACCTCATGACCGGGAACCACAGGCATCTTCACCATCGGGTTACCGCCACGGAAGGTATTCAGGTCACTGCCGCAAAAACC'


def test_load_fasta2(nucseq):
    seq_id, seq = next(fasta.load_fasta(nucseq))
    assert seq_id == u'contig-1467318'


def test_load_fasta3(nucseq):
    count = sum(1 for x in fasta.load_fasta(nucseq))
    assert count == 115


def test_load_fasta4(oneseq):
    count = sum(1 for x in fasta.load_fasta(oneseq))
    assert count == 1


def test_load_fasta5(oneseq):
    seq_id, seq = next(fasta.load_fasta(oneseq))
    assert len(seq) == 630


def test_write_fasta_sequence1(nucseq, tmpdir):

    seq_id, seq = next(fasta.load_fasta(nucseq))
    file_name = (tmpdir / 'test.fa').strpath
    file_handle = open_file(file_name, 'w')

    fasta.write_fasta_sequence(file_handle, seq_id, seq)
    file_handle.close()

    seq_idw, seqw = next(fasta.load_fasta(file_name))

    assert (seq_id, seq) == (seq_idw, seqw)


def test_write_fasta_sequence2(nucseq, tmpdir):

    file_name = (tmpdir / 'test.fa').strpath
    file_handle = open_file(file_name, 'w')

    for seq_id, seq in fasta.load_fasta(nucseq):
        fasta.write_fasta_sequence(file_handle, seq_id, seq)
    file_handle.close()

    count1 = sum(1 for x in fasta.load_fasta(nucseq))
    count2 = sum(1 for x in fasta.load_fasta(file_name))

    assert count1 == count2


def test_split_fasta_file1(nucseq, tmpdir):

    fasta.split_fasta_file(nucseq, (tmpdir / 'test{}.fa').strpath, 3)

    files = list(pathlib.Path(tmpdir.strpath).glob('*.fa'))

    assert len(files) == 3


def test_split_fasta_file2(nucseq, tmpdir):

    fasta.split_fasta_file(nucseq, (tmpdir / 'test{}.fa').strpath, 3)

    files = list(
        str(path)
        for path in pathlib.Path(tmpdir.strpath).glob('*.fa')
    )

    count1 = sum(1 for x in fasta.load_fasta(nucseq))
    count2 = sum(1 for x in fasta.load_fasta_files(files))

    assert count1 == count2


def test_load_fasta_rename1(oneseq):
    seq_id, seq = next(fasta.load_fasta_rename(oneseq))
    assert seq_id == 'M95099.1'


def test_load_fasta_rename2(oneseq):
    nam_func = lambda x: x.split('.')[0]
    seq_id, seq = next(fasta.load_fasta_rename(oneseq, name_func=nam_func))
    assert seq_id == 'M95099'
