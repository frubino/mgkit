import os.path
from nose import SkipTest

from mgkit.taxon import UniprotTaxonomy

base_dir = os.path.dirname(os.path.abspath(__file__))

data_dir = 'mg_data'
data_file = 'taxonomy.pickle'

try:
    TAXONOMY = UniprotTaxonomy(
        os.path.join(
            base_dir,
            data_dir,
            data_file
        )
    )
except IOError:
    TAXONOMY = None


def setup_taxon_data():
    if TAXONOMY is None:
        raise SkipTest(
            'No taxonomy data found: expecting file "{0}" in dir {1}'.format(
                data_file,
                data_dir
            )
        )
