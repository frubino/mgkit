"""
Metagenomics Framework
"""
from __future__ import print_function
import sys
import logging

from . import logger

__VERSION__ = "0.3.4"

__version__ = __VERSION__

LOG = logging.getLogger(__name__)

DEBUG = False
"Debug switch for a few functions"

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

PKG_NAME = 'MGKit'


def cite(file_handle=sys.stderr):
    """
    Print citation to the specified stream
    """
    print(
        LOGO,
        'MGKit Version: {0}'.format(__VERSION__),
        CITE,
        sep='\n',
        file=file_handle
    )


def check_version(version):
    if __version__ != version:
        LOG.warning(
            "This was tested with %s version %s (%s was found)",
            PKG_NAME,
            version,
            __version__
        )


class DependencyError(Exception):
    "Raised if an optional requirement is needed"
    def __init__(self, package):
        super(DependencyError, self).__init__(
            "The '{}' package is required".format(package)
        )
