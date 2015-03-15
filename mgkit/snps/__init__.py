"""
SNPs data package
"""

#these import are for "legacy" scripts
#from .funcs import combine_snps_in_dataframe_test
#combine_snps_in_dataframe = combine_snps_in_dataframe_test

#imported to support old pickled data and simpler access
from .classes import GeneSyn

#imported for convenience
from .conv_func import get_rank_dataframe, get_gene_map_dataframe
from .funcs import combine_sample_snps

import filter
import funcs
import conv_func
