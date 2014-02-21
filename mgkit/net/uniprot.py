"""
Contains function and constants for Uniprot access
"""

import urllib
import mgkit
import logging
from . import url_read

UNIPROT_MAP = 'http://www.uniprot.org/mapping/'
"URL to Uniprot mapping REST API"
UNIPROT_GET = 'http://www.uniprot.org/uniprot/'
"URL to Uniprot REST API"

LOG = logging.getLogger(__name__)


def get_sequences_by_ko(ko_id, taxonomy, contact=None, reviewed=True):
    """
    Gets sequences from Uniprot, restricting to the taxon id passed.

    :param str ko_id: KO id of the sequences to download
    :param int taxonomy: id of the taxon
    :param str contact: email address to be passed in the query (requested by
        Uniprot API)
    :param bool reviewed: if the sequences requested must be reviewed

    :return: string with the fasta file downloaded
    """
    params = urllib.urlencode(
        {
            'query': 'database:(type:ko {0}) AND taxonomy:{1}{2}'.format(
                ko_id, taxonomy, ' reviewed:yes' if reviewed else ''),
            'format': 'fasta',
            'limit': 200,
            'sort': 'score'
        }
    )
    if mgkit.DEBUG:
        LOG.debug("query: %s?%s", UNIPROT_GET, params)
        LOG.debug("request length %d", len(params))

    fasta = url_read(UNIPROT_GET, data=params, agent=contact)

    return fasta


def get_mappings(entry_ids, db_from='ID', db_to='EMBL', out_format='tab',
                 contact=None):
    """
    Gets mapping of genes using Uniprot REST API. The db_from and db_to values
    are the ones accepted by Uniprot API. The same applies to out_format, the
    only processed formats are 'list', which returns a list of the mappings
    (should be used with one gene only) and 'tab', which returns a dictionary
    with the mapping. All other values returns a string with the newline
    stripped.

    :param iterable entry_ids: iterable of ids to be mapped (there's a limit)
        to the maximum length of a HTTP request, so it should be less than 50
    :param str db_from: string that identify the DB for elements in entry_ids
    :param str db_to: string that identify the DB to which map entry_ids
    :param str out_format: format of the mapping; 'list' and 'tab' are processed
    :param str contact: email address to be passed in the query (requested
        Uniprot API)

    :return: tuple, dict or str depending on out_format value
    """

    if isinstance(entry_ids, str):
        entry_ids = [entry_ids]

    data = urllib.urlencode(
        {
            'from': db_from,
            'to': db_to,
            'query': ' '.join(entry_ids),
            'format': out_format
        }
    )

    mappings = url_read(UNIPROT_MAP, data=data, agent=contact)

    mappings = mappings.strip()

    if out_format == 'list':
        mappings = mappings.split('\n')
    elif out_format == 'tab':
        mapping_dict = {}
        mappings = mappings.split('\n')
        #delete first row 'From to'
        del mappings[0]

        for mapping in mappings:
            id_from, id_to = mapping.split('\t')
            if id_to == 'null':
                continue
            try:
                mapping_dict[id_from].append(id_to)
            except KeyError:
                mapping_dict[id_from] = [id_to]

        mappings = mapping_dict

    return mappings


def ko_to_mapping(ko_id, query, columns, contact=None):
    """
    Returns the mappings to the supplied KO. Can be used for any id, the
    query format is free as well as the columns returned. The only
    restriction is using a tab format, that is parsed.

    :param str ko_id: id used in the query
    :param str query: query passed to the Uniprot API, ko_id is replaced
        using :func:`str.format`
    :param str column: column used in the results table used to map the ids
    :param str contact: email address to be passed in the query (requested
        Uniprot API)

    .. note::

        each mapping in the column is separated by a ;

    """
    data = urllib.urlencode(
        {
            'query': query.format(ko_id),
            'format': 'tab',
            'columns': columns
        }
    )

    mappings = url_read(UNIPROT_GET, data=data, agent=contact)

    if mgkit.DEBUG:
        LOG.debug("query: %s?%s", UNIPROT_GET, data)
        LOG.debug("request length %d", len(data))

    mappings = mappings.split('\n')
    del mappings[0]

    categories = set()

    for map_line in mappings:

        mappings = [mapping.strip() for mapping in map_line.split(';')]
        if not mappings:
            continue

        categories.update(mappings)

    #in case an empty line is present
    try:
        categories.remove('')
    except KeyError:
        pass

    return categories
