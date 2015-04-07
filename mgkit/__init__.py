"""
Metagenomics Framework
"""
from __future__ import print_function
import sys

__VERSION__ = "0.2.0"

__version__ = __VERSION__

DEBUG = False
"Debug switch for a few functions"

from . import logger

LOGO = """
 _|      _|    _|_|_|  _|    _|  _|    _|
 _|_|  _|_|  _|        _|  _|        _|_|_|_|
 _|  _|  _|  _|  _|_|  _|_|      _|    _|
 _|      _|  _|    _|  _|  _|    _|    _|
 _|      _|    _|_|_|  _|    _|  _|      _|_|

"""

CITE = """
Rubino, F. and Creevey, C.J. (2014).
MGkit: Metagenomic Framework For The Study Of Microbial Communities.

Available at: http://figshare.com/articles/MGkit_Metagenomic_Framework_For_The_Study_Of_Microbial_Communities/1269288

[doi:10.6084/m9.figshare.1269288]
"""


def cite(file=sys.stderr):
    """
    Print citation to the specified stream
    """
    print(
        LOGO,
        'MGKit Version: {0}'.format(__VERSION__),
        CITE,
        sep='\n',
        file=file
    )
