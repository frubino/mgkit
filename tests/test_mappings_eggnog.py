import pytest
from conftest import skip_no_connection
import random
from mgkit.mappings.eggnog import NOGInfo

pytestmark = skip_no_connection


@pytest.fixture(scope='module')
def nog_info(eggnog_v3):
    eg = NOGInfo(**eggnog_v3)
    return eg


def test_load_members(eggnog_v3):
    eg = NOGInfo()
    eg.load_members(eggnog_v3['members'])
    assert eg._map_gene_nog
    assert eg._map_nog_gene


def test_load_description(eggnog_v3):
    eg = NOGInfo()
    eg.load_description(eggnog_v3['description'])
    assert eg._map_nog_desc


def test_load_funccat(eggnog_v3):
    eg = NOGInfo()
    eg.load_funccat(eggnog_v3['funccat'])
    assert eg._map_nog_func


def test_noginfo_get_nog_funccat1(nog_info):
    nog_id = 'NOG19649'
    assert nog_info.get_nog_funccat(nog_id)


def test_noginfo_get_nog_funccat2(nog_info):
    assert nog_info.get_nog_funccat('XXXXXXX') == set()


def test_noginfo_get_nogs_funccat1(nog_info):
    assert nog_info.get_nogs_funccat(
        ['NOG19649', 'NOG134552', 'NOG74662']
    ) == set(['S'])


def test_noginfo_get_nog_gencat(nog_info):
    nog_id = 'NOG74662'
    assert nog_info.get_nog_gencat(nog_id) == {'POORLY CHARACTERIZED'}


def test_noginfo_gene_nog1(nog_info):
    assert nog_info.get_gene_nog('NOG19649') is None


def test_noginfo_gene_nog2(nog_info):
    assert nog_info.get_gene_nog('5062.CADAORAP00008010') == {'NOG129254'}


def test_noginfo_get_gene_funccat1(nog_info):
    assert nog_info.get_gene_funccat('5062.CADAORAP00008010') == set()


def test_noginfo_get_gene_funccat2(nog_info):
    assert nog_info.get_gene_funccat('177437.HRM2_09030') == {'K'}


def test_noginfo_get_gene_funccat3(nog_info):
    assert nog_info.get_gene_funccat('XXXXXX') is None
