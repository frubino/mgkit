"Access EMBL Services"

import logging
import re
from . import url_read
import mgkit


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

EMBL_DBID = 'embl_cds'
"Default database id"

NOT_FOUND = r"Entry: .+? not found.\n"
"""
Appears in the resulting fasta (not tried on other formats) in the case that at
least one entry wasn't found.
"""

SUPPRESSED = "suppressed at the submitter's request on"

NONE_FOUND = r"ERROR 12.+?.\n?"
"Regular expression to check if no entry was found, used by :exc:`NoEntryFound`"

LOG = logging.getLogger(__name__)
"Log instance for this module"


def get_sequences_by_ids(embl_ids, contact=None, out_format='fasta', num_req=10,
                         embl_db=EMBL_DBID, compress=True, strict=False):
    """
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
        compress (bool): if True, the function tries to obtain a compressed
            response and decompress it on the fly
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

        seqs = url_read(request_url, compress=compress, agent=contact)

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

    #check for suppressed entries
    entries = '\n'.join(
        line for line in entries.split('\n')
        if line.find(SUPPRESSED) == -1
    )

    if error_re.search(entries):
        #if there's no sequences or we want all sequences to be downloaded
        if entries.find('>') == -1:
            raise NoEntryFound("No entry found")
        elif strict:
            raise EntryNotFound("At least one Entry was not found")
        #in case there's at least one sequence that was retrieved
        else:
            entries = ''.join(error_re.split(entries))
            LOG.debug("At least one Entry was not found")

    return entries


def dbfetch(embl_ids, db='embl', contact=None, out_format='seqxml', num_req=10):
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
