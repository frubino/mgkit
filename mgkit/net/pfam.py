"""
.. versionadded:: 0.2.3

This module defines routine to access Pfam information using a
network connection

"""
from . import url_open

PFAM_URL = "http://pfam.xfam.org/"


def get_pfam_families(key='id'):
    """
    .. versionadded:: 0.2.3

    Gets a dictionary with the accession/id/description of Pfam families
    from Pfam. This list can be accessed using the URL:
    http://pfam.xfam.org/families?output=text

    The output is a tab separated file where the fields are:

    * ACCESSION
    * ID
    * DESCRIPTION

    Arguments:
        key (str): if the value is *id*, the key of the dictionary is the ID,
            otherwise ID swaps position with ACCESSION (the new key)

    Returns:
        dict: by default the function returns a dictionary that uses the ID
        as key, while the value is a tuple (ACCESSION, DESCRIPTION). ID is the
        default because the :ref:`hmmer2gff` script output uses ID as *gene_id*
        value when using the HMM provided by Pfam
    """
    families = {}
    for line in url_open(PFAM_URL + "families?output=text", stream=True):
        line = line.decode('utf8').strip()
        if line.startswith('#') or (not line):
            continue
        acc, p_id, description = line.strip().split('\t')

        if key != 'id':
            p_id, acc = acc, p_id

        families[p_id] = (acc, description)

    return families
