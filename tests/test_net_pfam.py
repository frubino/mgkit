import pytest
from conftest import skipif
from mgkit.net.pfam import get_pfam_families

@skipif
def test_pfam_families_id():
    assert 'Protamine_P1' in get_pfam_families(key='id')


@skipif
def test_pfam_families_ac():
    assert 'PF07938' in get_pfam_families(key='ac')
