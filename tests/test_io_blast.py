import pytest

from mgkit.io import fasta
from mgkit.io.blast import parse_blast_tab, parse_uniprot_blast, \
    parse_accession_taxa_table

@pytest.fixture
def blast_tab_file(shared_datadir):
    return str(shared_datadir / 'blast-outfmt6.tsv')


@pytest.fixture
def blast_seq(shared_datadir):
    return dict(
        fasta.load_fasta_rename(str(shared_datadir / 'blast-seq.fa'))
    )


@pytest.fixture
def ncbi_taxa_file(shared_datadir):
    return str(shared_datadir / 'ncbi-taxa.tab')


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


def test_parse_blast_tab5(blast_tab_file):
    seq_id, hit = next(parse_blast_tab(
        blast_tab_file,
        key_func=lambda x: x.split('.')[0],
        ret_col=None)
    )
    assert len(hit) == 12


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


def test_parse_accession_taxa_table1(ncbi_taxa_file):
    assert next(parse_accession_taxa_table(ncbi_taxa_file, key=0, value=1,
                                           no_zero=False))[0] == 'KH113978'


def test_parse_accession_taxa_table2(ncbi_taxa_file):
    assert next(parse_accession_taxa_table(ncbi_taxa_file, key=0, value=1,
                                           no_zero=False))[1] == 0


def test_parse_accession_taxa_table3(ncbi_taxa_file):
    assert next(parse_accession_taxa_table(ncbi_taxa_file, key=0, value=1,
                                           no_zero=True))[0] == 'XM_006695435'


def test_parse_accession_taxa_table4(ncbi_taxa_file):
    acc_ids = ['XM_006695435']
    iterator = parse_accession_taxa_table(ncbi_taxa_file, key=0, value=1,
                                          acc_ids=acc_ids)
    iterator = list(acc_id for acc_id, taxon_id in iterator)
    assert acc_ids == iterator


def test_parse_accession_taxa_table5(ncbi_taxa_file):
    acc_ids = {'XM_006695435'}
    iterator = parse_accession_taxa_table(ncbi_taxa_file, key=0, value=1,
                                          acc_ids=acc_ids)
    iterator = list(acc_id for acc_id, taxon_id in iterator)
    assert list(acc_ids) == iterator
