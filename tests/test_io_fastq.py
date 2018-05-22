import pytest
from mgkit.io import fastq
from mgkit.io import open_file

@pytest.fixture
def fastq_file(shared_datadir):
    return str(shared_datadir / 'SP1.fq')


def test_load_fastq1(fastq_file):

    assert sum(1 for record in fastq.load_fastq(fastq_file)) == 250


def test_load_fastq2(fastq_file):

    header, seq, qual = next(fastq.load_fastq(fastq_file))

    assert seq == 'TTTCCGGGGCACATAATCTTCAGCCGGGCGC'


def test_load_fastq3(fastq_file):

    header, seq, qual = next(fastq.load_fastq(fastq_file))

    assert header == 'cluster_2:UMI_ATTCCG'


def test_load_fastq4(fastq_file):

    header, seq, qual = next(fastq.load_fastq(fastq_file))

    assert qual == '9C;=;=<9@4868>9:67AA<9>65<=>591'


def test_load_fastq4(fastq_file):

    header, seq, qual = next(fastq.load_fastq(fastq_file, num_qual=True))

    assert list(qual) == [24, 34, 26, 28, 26, 28, 27, 24, 31, 19, 23, 21, 23, 29, 24, 25, 21, 22, 32, 32, 27, 24, 29, 21, 20, 27, 28, 29, 20, 24, 16]


def test_write_fastq1(fastq_file, tmpdir):

    header, seq, qual = next(fastq.load_fastq(fastq_file))

    file_name = (tmpdir / 'test.fq').strpath

    file_handle = open_file(file_name, 'w')

    fastq.write_fastq_sequence(file_handle, header, seq, qual)
    file_handle.close()

    headerw, seqw, qualw = next(fastq.load_fastq(file_name))

    assert (header, seq, qual) == (headerw, seqw, qualw)


def test_write_fastq2(fastq_file, tmpdir):

    header, seq, qual = next(fastq.load_fastq(fastq_file, num_qual=True))

    file_name = (tmpdir / 'test.fq').strpath

    file_handle = open_file(file_name, 'w')

    fastq.write_fastq_sequence(file_handle, header, seq, qual)
    file_handle.close()

    headerw, seqw, qualw = next(fastq.load_fastq(file_name, num_qual=True))

    assert (header, seq, list(qual)) == (headerw, seqw, list(qualw))
