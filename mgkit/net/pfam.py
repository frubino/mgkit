"""
"""
from . import url_open

PFAM_URL = "http://pfam.xfam.org/"


def get_pfam_families(key='id'):
    """
    .. versionadded:: 0.2.3

    ACCESSION ID DESCRIPTION
    """
    families = {}
    for line in url_open(PFAM_URL + "families?output=text"):
        line = line.strip()
        if line.startswith('#') or (not line):
            continue
        acc, p_id, description = line.strip().split('\t')

        if key != 'id':
            p_id, acc = acc, p_id

        families[p_id] = (acc, description)

    return families
