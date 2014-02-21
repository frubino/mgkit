"""
Module containing constants for the filter package
"""

from ..snps import MIN_COV
from ..taxon import BLACK_LIST_IDS, is_ancestor

DEFAULT_SNP_FILTER = {
    'min_cov': MIN_COV,
    'black_list': BLACK_LIST_IDS,
    'func': is_ancestor
}
"""
Default filter options for filtering :class:`mgkit.snps.GeneSyn`
"""
