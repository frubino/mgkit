import pytest
from mgkit import kegg
from conftest import skip_no_connection

@pytest.fixture
def keggclient():
    return kegg.KeggClientRest()

pytestmark = skip_no_connection

def test_keggclient_link1(keggclient):
    assert keggclient.link_ids('rn', 'K00201') == {'K00201': ['R03015', 'R08060']}


def test_keggclient_list1(keggclient):
    assert 'cpd:C20660\tWybutosine in tRNA(Phe)\n' in keggclient.list_ids('cpd')


def test_keggclient_get1(keggclient):
    assert 'DBLINKS     PubChem: 172232382' in keggclient.get_entry('cpd:C20660')


def test_keggclient_get_names1(keggclient):
    assert 'M00002' in keggclient.get_ids_names('module')


def test_keggclient_get_names2(keggclient):
    assert 'md:M00002' in keggclient.get_ids_names('module', strip=False)


def test_cache_io(keggclient, tmpdir):
    query = keggclient.link_ids('rn', 'K00201')['K00201']
    file_name = tmpdir.join('cache.pickle').strpath
    keggclient.write_cache(file_name)

    cached = kegg.KeggClientRest(file_name)
    assert cached.cache == keggclient.cache


def test_cache_link_ids1(keggclient):
    query = keggclient.link_ids('rn', 'K00201')['K00201']
    assert query == keggclient.cache['link_ids']['rn']['K00201']


def test_find(keggclient):
    assert keggclient.find('CH4', 'compound') == {'C01438': 'Methane; CH4'}


def test_conv(keggclient):
    assert 'b0217' in keggclient.conv('ncbi-geneid', 'eco')


def test_get_ids_names(keggclient):
    assert 'K00001' in keggclient.get_ids_names('ko')
