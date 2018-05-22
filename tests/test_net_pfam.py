import pytest
from conftest import skip_no_connection
from mgkit.net.pfam import get_pfam_families

@skip_no_connection
def test_pfam_families_id():
    assert 'Protamine_P1' in get_pfam_families(key='id')


@skip_no_connection
def test_pfam_families_ac():
    assert 'PF07938' in get_pfam_families(key='ac')
