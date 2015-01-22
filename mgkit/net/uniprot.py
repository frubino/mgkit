"""
Contains function and constants for Uniprot access
"""

from __future__ import division

import urllib
import mgkit
import logging
import itertools
from . import url_read

UNIPROT_MAP = 'http://www.uniprot.org/mapping/'
"URL to Uniprot mapping REST API"

UNIPROT_GET = 'http://www.uniprot.org/uniprot/'
"URL to Uniprot REST API"

UNIPROT_TAXONOMY = 'http://www.uniprot.org/taxonomy/'
"URL to Uniprot REST API - Taxonomy"

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


def get_gene_info(gene_ids, columns, max_req=50, contact=None):
    """
    .. versionadded:: 0.1.12

    Get informations about a list of genes. it uses :func:`query_uniprot` to
    send the request and format the response in a dictionary.

    Arguments:
        gene_ids (iterable, str): gene id(s) to get informations for
        columns (list): list of columns
        max_req (int): number of maximum *gene_ids* per request
        contact (str): email address to be passed in the query (requested
            Uniprot API)

    Returns:
        dict: dictionary where the keys are the *gene_ids* requested and the
        values are dictionaries with the names of the *columns* requested as
        keys and the corresponding values, which can be lists if the values are
        are semicolon separated strings.

    Example:
        To get the taxonomy ids for some genes:

        >>> uniprot.get_gene_info(['Q09575', 'Q8DQI6'], ['organism-id'])
        {'Q09575': {'organism-id': '6239'}, 'Q8DQI6': {'organism-id': '171101'}}

    """
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]
    elif isinstance(gene_ids, set):
        gene_ids = list(gene_ids)

    if isinstance(columns, str):
        columns = [columns]

    infos = {}

    for index in range(0, len(gene_ids), max_req):

        LOG.info(
            "Querying uniprot ids (%d/%d)",
            index + max_req,
            len(gene_ids)
        )

        info_lines = query_uniprot(
            ' OR '.join(gene_ids[index:index+max_req]),
            columns=['id'] + columns,
            contact=contact
        )

        info_lines = info_lines.split('\n')

        del info_lines[0]

        for info_line in info_lines:
            info_line = info_line.strip()
            if not info_line:
                continue
            values = info_line.split('\t')

            gene_id = values[0]

            infos[gene_id] = dict(
                (
                    column,
                    value if (not value.endswith(';')) and (not value.endswith('; ')) and (not '; ' in value)
                    else [x.strip() for x in value.split(';') if x.strip()]
                )
                for column, value in zip(columns, values[1:])
            )

    return infos


def query_uniprot(query, columns=None, format='tab', limit=None, contact=None,
                  baseurl=UNIPROT_GET):
    """
    .. versionadded:: 0.1.12

    .. versionchanged:: 0.1.13
        added *baseurl* and made *columns* a default argument

    Queries Uniprot, returning the raw response in tbe format specified. More
    informations at the `page <http://www.uniprot.org/faq/28>`_

    Arguments:
        query (str): query to submit, as put in the input box
        columns (None, iterable): list of columns to return
        format (str): response format
        limit (int, None): number of entries to return or *None* to request all
            entries
        contact (str): email address to be passed in the query (requested
            Uniprot API)
        baseurl (str): base url for the REST API, can be either
            :data:`UNIPROT_GET` or :data:`UNIPROT_TAXONOMY`
    Returns:
        str: raw response from the query

    Example:
        To get the taxonomy ids for some genes:

        >>> uniprot.query_uniprot('Q09575 OR Q8DQI6', ['id', 'organism-id'])
        'Entry\\tOrganism ID\\nQ8DQI6\\t171101\\nQ09575\\t6239\\n'

    .. warning::

        because of limits in the length of URLs, it's advised to limit the
        length of the query string.

    """
    data = {
        'query': query,
        'format': format
    }

    if limit is not None:
        data['limit'] = limit

    data = urllib.urlencode(data)

    if columns is not None:
        data += "&columns={0}".format(
            ','.join(
                urllib.quote(column)
                for column in columns
            )
        )

    if mgkit.DEBUG:
        LOG.debug("query: %s?%s", baseurl, data)
        LOG.debug("request length %d", len(data))

    return url_read(baseurl, data, agent=contact)


def parse_uniprot_response(data, simple=True):
    """
    .. versionadded:: 0.1.12

    Parses raw response from a Uniprot query (tab format only) from functions
    like :func:`query_uniprot` into a dictionary. It requires that the first
    column is the entry id (or any other unique id).

    Arguments:
        data (str): string response from Uniprot
        simple (bool): if True and the number of columns is 1, the dictionary
            returned has a simplified structure

    Returns:
        dict: The format of the resulting dictionary is
        entry_id -> {column1 -> value, column2 -> value, ..} unless there's
        only one column and *simple* is True, in which case the value is
        equal to the value of the only column.
    """
    data = data.splitlines()

    columns = [x.lower() for x in data[0].split('\t')[1:]]

    del data[0]

    parsed_data = {}

    for line in data:
        line = line.split('\t')
        entry_id = line[0]

        if (len(columns) == 1) and simple:
            parsed_data[entry_id] = line[1]
        else:
            parsed_data[entry_id] = dict(
                itertools.izip(columns, line[1:])
            )

    return parsed_data


def get_ko_to_eggnog_mappings(ko_ids, contact=None):
    """
    .. versionadded:: 0.1.14

    It's not possible to map in one go KO IDs to eggNOG IDs via the API in
    Uniprot. This function uses :func:`query_uniprot` to get all Uniprot IDs
    requested and the return a dictionary with all their eggNOG IDs they map to.

    Arguments:
        ko_ids (iterable): an iterable of KO IDs
        contact (str): email address to be passed in the query (requested
            Uniprot API)

    Returns:
        dict: The format of the resulting dictionary is
        ko_id -> {eggnog_id1, ..}
    """

    data = query_uniprot(
        "database:(type:ko AND ({}))".format(' OR '.join(list(ko_ids))),
        columns=['database(KO)', 'database(EGGNOG)'],
        contact=contact
    )

    data = data.splitlines()

    del data[0]

    parsed_data = {}

    for line in data:
        ko_ids, eggnog_ids = line.split('\t')
        ko_ids = ko_ids.split(';')

        for ko_id in ko_ids:
            if not ko_id:
                continue
            if not eggnog_ids:
                continue

            for eggnog_id in eggnog_ids.split(';'):
                if not eggnog_id:
                    continue

                try:
                    parsed_data[ko_id].add(eggnog_id)
                except KeyError:
                    parsed_data[ko_id] = set([eggnog_id])

    return parsed_data


def get_uniprot_ec_mappings(gene_ids, contact=None):
    """
    .. versionadded:: 0.1.14

    Shortcut to download EC mapping of Uniprot IDs. Uses :func:`get_gene_info`
    passing the correct column (*ec*).

    """
    return get_gene_info(
        gene_ids,
        columns=['ec'],
        contact=contact,
        max_req=100
    )
