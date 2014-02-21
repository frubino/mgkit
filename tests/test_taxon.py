from nose.tools import *
import random

from mgkit.taxon import is_ancestor

import taxon_data


@with_setup(setup=taxon_data.setup_taxon_data)
def test_is_ancestor1():
    #prevotella ruminicola, 839
    taxon_id = 839
    anc_id = taxon_data.TAXONOMY.find_by_name('prevotella')[0]
    eq_(
        is_ancestor(taxon_data.TAXONOMY, taxon_id, anc_id),
        True
    )


@with_setup(setup=taxon_data.setup_taxon_data)
def test_is_ancestor2():
    #any prevotella species
    taxon_id = random.choice(
        [
            taxon_obj.taxon_id
            for taxon_obj in taxon_data.TAXONOMY
            if 'prevotella ' in taxon_obj.s_name
        ]
    )
    anc_id = taxon_data.TAXONOMY.find_by_name('clostridium')[0]
    eq_(
        is_ancestor(taxon_data.TAXONOMY, taxon_id, anc_id),
        False
    )
