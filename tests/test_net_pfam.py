import pytest

from mgkit.net.pfam import get_pfam_families

def test_pfam_families_id():
    assert 'Protamine_P1' in get_pfam_families(key='id')


def test_pfam_families_ac():
    assert 'PF07938' in get_pfam_families(key='ac')
