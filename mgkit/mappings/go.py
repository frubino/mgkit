
"""
Module containing classes and functions to deal with Gene Ontology data
"""
from future.utils import viewitems
from requests.exceptions import HTTPError
import logging
from .. import kegg
from ..utils import dictionary as dict_utils

LOG = logging.getLogger(__name__)
