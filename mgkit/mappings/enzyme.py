"""
.. versionadded:: 0.1.14

EC mappings
"""

import re
from ..io import open_file

LEVEL1_NAMES = {
    1: 'oxidoreductases',
    2: 'transferases',
    3: 'hydrolases',
    4: 'lyases',
    5: 'isomerases',
    6: 'ligases'
}
"""
Top level classification names
"""

ENZCLASS_REGEX = r"^(\d)\. ?([\d-]+)\. ?([\d-]+)\. ?([\d-]+) +(.+)\."
"""
Used to get the description for the higher level enzyme classes from the file
*enzclass.txt* on `expasy <http://expasy.org>`_
"""


def parse_expasy_file(file_name):
    """
    Used to load enzyme descriptions from the file *enzclass.txt* on
    `expasy <http://expasy.org>`_
    """
    labels = {}

    for line in open_file(file_name, mode='r'):
        match = re.search(ENZCLASS_REGEX, line)

        if match is None:
            continue

        ec = '.'.join(value for value in match.groups()[:4] if value != '-')

        labels[ec] = match.group(5)

    return labels


def get_enzyme_level(ec, level=4):
    """
    .. versionadded:: 0.1.14

    Returns an enzyme class at a specific level , between 1 and 4 (by default
    the most specific, 4)

    Arguments:
        ec (str): a string representing an EC number (e.g. 1.2.4.10)
        level (int): from 1 to 4, to get a different level specificity of in
            the enzyme classification

    Returns:
        str: the EC number at the requested specificity

    Example:
        >>> from mgkit.mappings.enzyme import get_enzyme_level
        >>> get_enzyme_level('1.1.3.4', 1)
        '1'
        >>> get_enzyme_level('1.1.3.4', 2)
        '1.1'
        >>> get_enzyme_level('1.1.3.4', 3)
        '1.1.3'
        >>> get_enzyme_level('1.1.3.4', 4)
        '1.1.3.4'
    """
    return '.'.join(ec.split('.')[:level])


def change_mapping_level(ec_map, level=3):
    """
    .. versionadded:: 0.1.14

    Given a dictionary, whose values are dictionaries, in which a key is named
    *ec* and its value is an iterable of EC numbers, returns an iterator that
    can be used to build a dictionary with the same top level keys and the
    values are sets of the transformed EC numbers.

    Arguments:
        ec_map (dict): dictionary generated by
            :func:`mgkit.net.uniprot.get_gene_info`
        level (int): number from 1 to 4, to specify the level of the mapping,
            passed to :func:`get_enzyme_level`

    Yields:
        tuple: a tuple (gene_id, set(ECs)), which can be passed to *dict* to
        make a dictionary

    Example:
        >>> from mgkit.net.uniprot import get_gene_info
        >>> from mgkit.mappings.enzyme import change_mapping_level
        >>> ec_map = get_gene_info('Q9HFQ1', columns='ec')
        {'Q9HFQ1': {'ec': '1.1.3.4'}}
        >>> dict(change_mapping_level(ec_map, level=2))
        {'Q9HFQ1': {'1.1'}}

    """
    for gene_id, ecdict in ec_map.iteritems():
        try:
            ec_list = ecdict['ec']
        except KeyError:
            continue

        if isinstance(ec_list, str):
            ec_list = [ec_list]

        yield gene_id, set(
            get_enzyme_level(ec_id, level=level) for ec_id in ec_list
        )
