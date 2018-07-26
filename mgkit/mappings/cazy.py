"""
Module containing classes and functions to deal with CaZy data
"""

from requests.exceptions import HTTPError
import logging
from .. import kegg

CAZY_FAMILIES = {
    'GH': 'Glycoside Hydrolase',
    'GT': 'GlycosylTransferase',
    'PL': 'Polysaccharide Lyase',
    'CE': 'Carbohydrate Esterase',
    'CBM': 'Carbohydrate-Binding Module'
}
"CaZy families"

LOG = logging.getLogger(__name__)
