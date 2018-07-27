import pytest
from mgkit.net.uniprot import get_sequences_by_ko, get_mappings, get_gene_info, \
    query_uniprot
from conftest import skip_no_connection


@skip_no_connection
def test_get_sequences_by_ko():
    assert get_sequences_by_ko('K00001', 2)[0] == '>'


@skip_no_connection
def test_get_mappings1():
    assert 'A0A023P6Q7' in \
        get_mappings('K00001', db_from='KO_ID', db_to='ACC')['K00001']


def test_get_mappings2():
    assert 'A0A023P6Q7' in \
        get_mappings('K00001', db_from='KO_ID', db_to='ACC', out_format='list')


# @pytest.mark.skip(reason='There are discrepancies in Uniprot return values')
@skip_no_connection
def test_get_gene_info():
    assert get_gene_info(['Q09575', 'Q8DQI6'], ['organism-id']) == \
        {'Q09575': {'organism-id': '6239'}, 'Q8DQI6': {'organism-id': '171101'}}


@skip_no_connection
def test_query_uniprot():
    result = query_uniprot('Q09575 OR Q8DQI6', ['id', 'organism-id'])
    assert ('Entry' in result) and ('Organism ID' in result)
