"Access EMBL Services"

from builtins import range
import logging
import re
import gzip
from . import url_read
import requests
import mgkit
try:
    from cStringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO


class NoEntryFound(Exception):
    """
    Raised if no sequences where found by :func:`get_sequences_by_ids`, the
    check is based on the :data:`NONE_FOUND` variable.
    """
    pass


class EntryNotFound(Exception):
    """
    Raised if at least one entry was not found by :func:`get_sequences_by_ids`.
    :data:`NOT_FOUND` is used to check if any entry wasn't downloaded.
    """
    pass


URL_REST = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/"
"Base URL for EMBL DBFetch REST API"

URL_DATAWAREHOUSE = "http://www.ebi.ac.uk/ena/data/warehouse/search?"

EMBL_DBID = 'embl_cds'
"Default database id"

NOT_FOUND = r"Entry: .+? not found.\n"
"""
Appears in the resulting fasta (not tried on other formats) in the case that at
least one entry wasn't found.
"""

SUPPRESSED = "suppressed at the submitter's request on"

NONE_FOUND = r"ERROR 12.+?.\n?"
"""
Regular expression to check if no entry was found, used by :exc:`NoEntryFound`
"""

LOG = logging.getLogger(__name__)
"Log instance for this module"


def get_sequences_by_ids(embl_ids, contact=None, out_format='fasta',
                         num_req=10, embl_db=EMBL_DBID, strict=False):
    """
    .. versionchanged:: 0.3.4
        removed *compress* as it's bases on the `requests` package

    Downloads entries using EBI REST API. It can download one entry at a
    time or accept an iterable and all sequences will be downloaded in batches
    of at most num_req.

    It's fairly general, so can be customised, from the DB used to the output
    format: all batches are simply concatenate.

    .. note::

        There are some checks on the some errors reported by the EMBL api, but
        not documented, in particular two errors, which are just reported as
        text lines in the fasta file (the only one tested at this time).

        The are two possible cases:

        * if no entry was found :exc:`NoEntryFound` will be raised.
        * if at least one entry wasn't found:

          * if strict is False (the default) the error will be just logged as a
            debug message
          * if strict is True :exc:`EntryNotFound` is raised

    Args:
        embl_ids (iterable, str): list of ids to download
        contact (str): email address to be passed in the query
        format (str): format of the entry
        num_req (int): number of entries to download with each request
        embl_db (str): db to which the ids refer to
        strict (bool): if True, a check on the number of entries retrieved is
            performed

    Returns:
        str: the entries requested

    Raises:
        EntryNotFound: if at least an entry was not found
        NoEntryFound: if NO entry were found

    .. warning::

        The number of sequences that can be downloaded at a time is 11, it
        seems, since the returned sequences for each request was at most 11. I
        didn't find any mention of this in the API docs, but it may be a
        restriction that's temporary.

    """
    if isinstance(embl_ids, str):
        embl_ids = [embl_ids]
    elif isinstance(embl_ids, (set, dict)):
        embl_ids = list(embl_ids)

    splitter = re.compile(NOT_FOUND)
    error_re = re.compile(NONE_FOUND)

    entries = []

    for idx in range(0, len(embl_ids), num_req):

        if len(embl_ids) > num_req:
            LOG.debug("Downloading ids - range %d-%d", idx + 1, idx + num_req)

        request_params = "{db_id}/{entry_id}/{format}".format(
            db_id=embl_db,
            entry_id=','.join(embl_ids[idx:idx+num_req]),
            format=out_format
        )

        request_url = URL_REST + request_params

        seqs = url_read(request_url, agent=contact)

        # it seems that they add at the end of each request two newlines
        # (should be just one)
        seqs = seqs.replace('\n\n', '\n')

        if splitter.search(seqs):
            if strict:
                raise EntryNotFound("At least one Entry was not found")
            LOG.debug("At least one Entry was not found")
            entries.extend(splitter.split(seqs))
        else:
            entries.append(seqs)

    entries = ''.join(entries)

    # check for suppressed entries
    entries = '\n'.join(
        line for line in entries.split('\n')
        if line.find(SUPPRESSED) == -1
    )

    if error_re.search(entries):
        # if there's no sequences or we want all sequences to be downloaded
        if entries.find('>') == -1:
            raise NoEntryFound("No entry found")
        elif strict:
            raise EntryNotFound("At least one Entry was not found")
        # in case there's at least one sequence that was retrieved
        else:
            entries = ''.join(error_re.split(entries))
            LOG.debug("At least one Entry was not found")

    return entries


def dbfetch(embl_ids, db='embl', contact=None, out_format='seqxml',
            num_req=10):
    """
    .. versionadded:: 0.1.12

    Function that allows to use dbfetch service (REST). More information on the
    output formats and the database available at the
    `service page <http://www.ebi.ac.uk/Tools/dbfetch/syntax.jsp>`_

    Arguments:
        embl_ids (str, iterable): list or single sequence id to retrieve
        db (str): database from which retrieve the sequence data
        contact (str): email contact to use as per EMBL guidlines
        out_format (str): output format, depends on database
        num_req (int): number of ids per request

    Returns:
        list: a list with the results from each request sent. Each request sent
        has a maximum number *num_req* of ids, so the number of items in the
        list depends by the number of ids in *embl_ids* and the value of
        *num_req*.
    """
    if isinstance(embl_ids, str):
        embl_ids = [embl_ids]
    elif isinstance(embl_ids, (set, dict)):
        embl_ids = list(embl_ids)

    entries = []

    for idx in range(0, len(embl_ids), num_req):

        if len(embl_ids) > num_req:
            LOG.info("Downloading ids - range %d-%d", idx + 1, idx + num_req)

        request_params = "{db_id}/{entry_id}/{format}?style=raw".format(
            db_id=db,
            entry_id=','.join(embl_ids[idx:idx+num_req]),
            format=out_format
        )

        request_url = URL_REST + request_params

        if mgkit.DEBUG:
            LOG.debug(request_url)

        entries.append(url_read(request_url, agent=contact))

    return entries


def datawarehouse_search(query, domain='sequence', result='sequence_release',
                         display='fasta', offset=0, length=100000,
                         contact=None, download='gzip', url=URL_DATAWAREHOUSE,
                         fields=None):
    """
    .. versionchanged:: 0.2.3
        added *fields* parameter to retrieve tab separated information

    .. versionadded:: 0.1.13

    Perform a datawarehouse search on EMBL dbs. Instructions on the query
    language used to query the datawarehouse are available at `this page
    <http://www.ebi.ac.uk/ena/about/browser#data_warehouse>`_ with more details
    about the databases domains `at this page
    <http://www.ebi.ac.uk/ena/data/warehouse/usage>`_

    Arguments:
        query (str): query for the search enging
        domain (str): database domain to search
        result (str): domain result requested
        display (str): display option (format to retrieve the entries)
        offset (int): the offset of the search results, defaults to the first
        length (int): number of results to retrieve at the specified offset
            and the limit is automatically set a 100,000 records for query
        contact (str): email of the user
        download (str): type of response. Gzip responses are automatically
            decompressed
        url (str): base URL for the resource
        fields (None, iterable): must be an iterable of fields to be returned
            if display is set to *report*

    Returns:
        str: the raw request

    Examples:
        Querying EMBL for all sequences of type rRNA of the  *Clostridium*
        genus. Only from the EMBL release database in fasta format:

        >>> query = 'tax_tree(1485) AND mol_type="rRNA"'
        >>> result = 'sequence_release'
        >>> display = 'fasta'
        >>> data = embl.datawarehouse_search(query, result=result,
        ... display=display)
        >>> len(data)
        35919

        Each entry taxon_id from the same data can be retrieved using *report*
        as the *display* option and *fields* an iterable of fields to just
        ('accession', tax_id'):

        >>> query = 'tax_tree(1485) AND mol_type="rRNA"'
        >>> result = 'sequence_release'
        >>> display = 'report'
        >>> fields = ('accession', 'tax_id')
        >>> data = embl.datawarehouse_search(query, result=result,
            display=display, fields=fields)

    """

    params = {
        'query': query,
        'domain': domain,
        'result': result,
        'display': display,
        'offset': offset,
        'length': length,
        'download': download
    }

    # if display='report', fields can be used to display, but download needs to
    # be 'txt'
    if display == 'report':
        params['download'] = download = 'txt'
        try:
            params['fields'] = ','.join(fields)
        except TypeError:
            raise ValueError(
                "the 'fields' parameter must be an iteratble if display is "
                "equal to 'report'"
                )

    stream = True if download == 'gzip' else False
    request = requests.get(url, params, headers={'user-agent': contact}, stream=stream)
    if download == 'gzip':
        data = StringIO(request.raw.read())
        data = gzip.GzipFile(mode='rb', fileobj=data).read()
    else:
        data = request.text

    return data
