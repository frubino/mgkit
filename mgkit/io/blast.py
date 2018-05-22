"""
Blast routines and parsers

"""
import sys
import logging
from builtins import range, zip
from future.utils import viewitems
from . import gff
from . import open_file
from ..utils.common import deprecated

NUM_LINES = 10 ** 6

LOG = logging.getLogger(__name__)


def add_blast_result_to_annotation(annotation, gi_taxa_dict, taxonomy,
                                   threshold=60):
    """
    Adds blast information to a GFF annotation.

    :param annotation: GFF annotation object
    :param dict gi_taxa_dict: dictionary returned by
        :func:`parse_gi_taxa_table`.
    :param taxonomy: Uniprot taxonomy, used to add the taxon name to the
        annotation
    """

    try:
        hit = gi_taxa_dict[annotation.attributes.ko_idx]
    except KeyError:
        return

    if hit.score <= threshold:
        return

    try:
        # skips if the ID is not in the Uniprot Taxonomy
        annotation.attributes.blast_taxon = taxonomy[hit.taxon_id].s_name
        annotation.attributes.blast_taxon_idx = hit.taxon_id
    except KeyError:
        return


def parse_blast_tab(file_handle, seq_id=0, ret_col=(0, 1, 2, 6, 7, 11),
                    key_func=None, value_funcs=None):
    """
    .. versionadded:: 0.1.12

    Parses blast output tab format, returning for each line a key (the query
    id) and the columns requested in a tuple.

    Arguments:
        file_handle (file): file name or file handle for the blast ouput
        seq_id (int): index for the column which has the query id
        ret_col (list, None): list of indexes for the columns to be returned or
            *None* if all columns must be returned
        key_func (None, func): function to transform the query id value in the
            key returned. If *None*, the query id is used
        value_funcs (None, list): list of functions to transform the value of
            all the requested columns. If *None* the values are not converted

    Yields:
        tuple: iterator of tuples with the first element being the query id
        after key_func is applied, if requested and the second element of
        the tuple is a tuple with the requested columns *ret_col*

    .. table:: BLAST+ used with `-outfmt 6`, default columns

        +--------------+--------------------------------+
        | column index |          description           |
        +==============+================================+
        |            0 | query name                     |
        +--------------+--------------------------------+
        |            1 | subject name                   |
        +--------------+--------------------------------+
        |            2 | percent identities             |
        +--------------+--------------------------------+
        |            3 | aligned length                 |
        +--------------+--------------------------------+
        |            4 | number of mismatched positions |
        +--------------+--------------------------------+
        |            5 | number of gap positions        |
        +--------------+--------------------------------+
        |            6 | query sequence start           |
        +--------------+--------------------------------+
        |            7 | query sequence end             |
        +--------------+--------------------------------+
        |            8 | subject sequence start         |
        +--------------+--------------------------------+
        |            9 | subject sequence end           |
        +--------------+--------------------------------+
        |           10 | e-value                        |
        +--------------+--------------------------------+
        |           11 | bit score                      |
        +--------------+--------------------------------+

    """

    if key_func is None:
        key_func = lambda x: x

    if value_funcs is None:
        value_funcs = [lambda x: x for y in range(len(ret_col))]

    if ret_col is None:
        ret_col = list(range(12))

    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading BLAST results from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    lineno = 0
    comments = 0

    for lineno, line in enumerate(file_handle):
        line = line.decode('ascii')
        if line.startswith('#'):
            comments += 1
            continue

        line = line.strip()

        if not line:
            continue

        line = line.split('\t')
        key = key_func(line[seq_id])
        values = tuple(
            func(line[index])
            for index, func in zip(ret_col, value_funcs)
        )

        yield key, values

    LOG.info('Read %d BLAST records', lineno + 1 - comments)


def parse_uniprot_blast(file_handle, bitscore=40, db='UNIPROT-SP', dbq=10,
                        name_func=None, feat_type='CDS', seq_lengths=None):
    """
    .. versionadded:: 0.1.12

    .. versionchanged:: 0.1.13
        added *name_func* argument

    .. versionchanged:: 0.2.1
        added *feat_type*

    .. versionchanged:: 0.2.3
        added *seq_lengths* and added subject *start* and *end* and *e-value*

    Parses BLAST results in tabular format using :func:`parse_blast_tab`,
    applying a basic bitscore filter. Returns the annotations associated with
    each BLAST hit.

    Arguments:
        file_handle (str, file): file name or open file handle
        bitscore (int, float): the minimum bitscore for an annotation to be
            accepted
        db (str): database used
        dbq (int): an index indicating the quality of the sequence database
            used; this value is used in the filtering of annotations
        name_func (func): function to convert the name of the database
            sequences. Defaults to `lambda x: x.split('|')[1]`, which can be
            be used with fasta files provided by Uniprot
        feat_type (str): feature type in the GFF
        seq_lengths (dict): dictionary with the sequences lengths, used to
            deduct the frame of the '-' strand

    Yields:
        Annotation: instances of :class:`mgkit.io.gff.Annotation` instance of
        each BLAST hit.
    """

    if name_func is None:
        name_func = lambda x: x.split('|')[1]

    ret_col = (0, 1, 2, 6, 7, 8, 9, 10, 11)

    # the second function extract the Uniprot ID from the sequence header
    value_funcs = (
        str,
        name_func,
        float,
        int,
        int,
        int,
        int,
        float,
        float
    )

    for seq_id, hit in parse_blast_tab(file_handle, ret_col=ret_col,
                                       value_funcs=value_funcs):
        if hit[-1] < bitscore:
            continue

        seq_len = None if seq_lengths is None else seq_lengths[seq_id]

        yield gff.from_nuc_blast(hit, db=db, dbq=dbq, feat_type=feat_type,
                                 seq_len=seq_len)


def parse_fragment_blast(file_handle, bitscore=40.0):
    """
    .. versionadded:: 0.1.13

    Parse the output of a BLAST output where the sequences are the single
    annotations, so the sequence names are the *uid* of the annotations.

    The only returned values are the best hits, maxed by bitscore and identity.

    Arguments:
        file_handle (str, file): file name or open file handle
        bitscore (float): minimum bitscore for accepting a hit

    Yields:
        tuple: a tuple whose first element is the *uid* (the sequence name) and
        the second is the a list of tuples whose first element is the GID (NCBI
        identifier), the second one is the identity and the third is the
        bitscore of the hit.

    """

    value_funcs = (lambda x: x.split('|')[1], float, float)

    iterator = parse_blast_tab(
        file_handle,
        ret_col=(1, 2, 11),
        value_funcs=value_funcs
    )

    uidmap = {}

    for uid, hit in iterator:

        if hit[-1] < bitscore:
            continue

        try:
            uidmap[uid].append(hit)
        except KeyError:
            uidmap[uid] = [hit]

    for uid, hits in viewitems(uidmap):
        # returns the hit with the max bitscore and max identity
        yield uid, hits


def parse_accession_taxa_table(file_handle, acc_ids=None, key=1, value=2,
                               num_lines=NUM_LINES, no_zero=True):
    """
    .. versionadded:: 0.2.5

    .. versionchanged:: 0.3.0
        added *no_zero*

    This function superseeds :func:`parse_gi_taxa_table`, since NCBI is
    deprecating the GIDs in favor of accessions like *X53318*. The new file can
    be found at the NCBI `ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid`,
    for DNA sequences (*nt* DB) *nucl_gb.accession2taxid.gz*.

    The file contains 4 columns, the first one is the accession without its
    version, the second one includes the version, the third column is the
    taxonomic identifier and the fourth is either the old GID or **na**.

    The column used as key is the *second*, since by default the fasta headers
    used in NCBI DBs use the versioned identifier. To use the GID as key, the
    *key* parameter can be set to 3, but if no identifier is found (*na* as per
    the file README), the line is skipped.

    Arguments:
        file_handle (str, file): file name or open file handle
        acc_ids (None, list): if it's not `None` only the keys included in the
            passed `acc_ids` list will be returned
        key (int): 0-based index for the column to use as accession. Defaults
            to the versioned accession that is used in GenBank fasta files.
        num_lines (None, int): number of which a message is logged. If None,
            no message is logged
        no_zero (bool): if True (default) a key with taxon_id of 0 is not yield

    .. note::

        GIDs are being phased out in September 2016:
        http://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/

    """

    if (sys.version_info[0] == 2) and isinstance(file_handle, unicode):
        file_handle = open_file(file_handle, 'rb')
    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading taxonomic information from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    if acc_ids is not None:
        acc_ids = set(acc_ids)

    zero_acc = 0

    for idx, line in enumerate(file_handle):
        line = line.decode('ascii')
        # skip header
        if line.startswith('accession'):
            continue

        if (num_lines is not None) and ((idx + 1) % num_lines == 0):
            LOG.info("Parsed %d lines", idx + 1)

        line = line.strip().split('\t')

        acc_id = line[key]

        if (acc_ids is not None) and (acc_id not in acc_ids):
            continue

        # this is the case if the column used as key is the fourth (GID) and
        # the GID is not available for that accession. The best case is to skip
        if acc_id.lower() == 'na':
            continue

        taxon_id = int(line[value])

        # in cases where 0 is the organism (doesn't exist in the taxonomy)
        if no_zero and (taxon_id == 0):
            # LOG.warning("accession '%s' has taxon_id 0", acc_id)
            zero_acc += 1
            continue

        yield acc_id, taxon_id

    if no_zero and (zero_acc > 0):
        LOG.warning("%d accessions have taxon_id 0", zero_acc)

    LOG.info("Parsed %d lines", idx + 1)
