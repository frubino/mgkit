"""
.. versionadded:: 0.1.14

EC mappings
"""

import re
from ..io import open_file

enzclass_regex = r"^(\d)\. ?([\d-]+)\. ?([\d-]+)\. ?([\d-]+) +(.+)\."
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
        match = re.search(enzclass_regex, line)

        if match is None:
            continue

        ec = '.'.join(value for value in match.groups()[:4] if value != '-')

        labels[ec] = match.group(5)

    return labels


def get_enzyme_level(ec, level=4):
    """
    Returns an enzyme class at a specific level , between 1 and 4 (by default
    the most specific, 4)
    """
    return '.'.join(ec.split('.')[:level])
