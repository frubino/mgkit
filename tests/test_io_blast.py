import pytest

from mgkit.io import fasta
from mgkit.io.blast import parse_blast_tab, parse_uniprot_blast

@pytest.fixture
def blast_tab_file(shared_datadir):
    return str(shared_datadir / 'blast-outfmt6.tsv')


@pytest.fixture
def blast_seq(shared_datadir):
    return dict(
        fasta.load_fasta_rename(str(shared_datadir / 'blast-seq.fa'))
    )


@pytest.fixture
def blast_lengths(blast_seq):
    return {
        seq_id: len(blast_seq[seq_id])
        for seq_id in blast_seq
    }

def test_parse_blast_tab1(blast_tab_file):
    seq_id, hit = next(parse_blast_tab(blast_tab_file))
    assert seq_id == 'M95099.1'


def test_parse_blast_tab2(blast_tab_file):
    seq_id, hit = next(parse_blast_tab(blast_tab_file))
    assert len(hit) == 6


def test_parse_blast_tab3(blast_tab_file):
    seq_id, hit = next(parse_blast_tab(blast_tab_file))
    assert hit[-1] ==  '14520'


def test_parse_blast_tab4(blast_tab_file):
    seq_id, hit = next(parse_blast_tab(blast_tab_file, key_func=lambda x: x.split('.')[0]))
    assert seq_id == 'M95099'


def test_parse_uniprot_blast1(blast_tab_file):
    key_func = lambda x: x.split('.')[0]

    a = next(parse_uniprot_blast(blast_tab_file, name_func=key_func))

    assert a.gene_id == 'M95099'


def test_parse_uniprot_blast2(blast_tab_file):
    key_func = lambda x: x.split('.')[0]

    a = next(parse_uniprot_blast(blast_tab_file, name_func=key_func))

    assert a.bitscore == 14520


def test_parse_uniprot_blast3(blast_tab_file):
    key_func = lambda x: x.split('.')[0]

    bts = [
        a.bitscore
        for a in parse_uniprot_blast(blast_tab_file, name_func=key_func, bitscore=170)
    ]

    assert bts == [14520]


def test_parse_uniprot_blast4(blast_tab_file, blast_lengths):
    key_func = lambda x: x.split('.')[0]

    a = next(parse_uniprot_blast(blast_tab_file, name_func=key_func, seq_lengths=blast_lengths))

    assert a.phase == 0
