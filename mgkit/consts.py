"""
Module containing constants for the filter package
"""

import mgkit.taxon

MIN_COV = 4
"Minumum coverage required in some functions."

MIN_NUM = 10
"Used to set the minimum number of replicates for some functions"

BLACK_LIST = [
    'bos',
    'pecora',
    'lolium',
    'streptophyta',
    'oryza',
    'fabales',
    'poaceae',
    'metazoa',
    'chlorophyta'
]
"""
Default taxa black list, includes all taxa names that are to be excluded from
some analysis.
"""

BLACK_LIST_IDS = [
    903,    # bos
    35500,  # pecora
    4520,   # lolium
    35493,  # streptophyta
    4527,   # oryza
    72025,  # fabales
    4479,   # poaceae
    33208,  # metazoa
    3041,   # chlorophyta
]

DEFAULT_SNP_FILTER = {
    'min_cov': MIN_COV,
    'black_list': BLACK_LIST_IDS,
    'include_only': [
        mgkit.taxon.ARCHAEA,
        mgkit.taxon.BACTERIA,
        mgkit.taxon.FUNGI,
    ] + mgkit.taxon.PROTISTS.values(),
    'func': mgkit.taxon.is_ancestor
}
"""
Default filter options for filtering :class:`mgkit.snps.GeneSyn`
"""
