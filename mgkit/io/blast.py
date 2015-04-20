"""
Blast routines and parsers

"""

import logging
import collections
from . import gff
from . import open_file
from ..utils.common import deprecated

NUM_LINES = 10 ** 6

LOG = logging.getLogger(__name__)


@deprecated
def _parse_blast_tab(f_handle, gid_col=1, score_col=11, num_lines=NUM_LINES):
    """
    .. deprecated:: 0.1.12
        Use the more general version :func:`parse_blast_tab`

    Parses blast output tab format, taking the best hit (first one) for each
    KO id (ko_idx). Returns a dictionary with ko_idx as key and the GID and bit
    score as value.

    :param file f_handle: file handle for the blast ouput
    :param int gid_col: index for the column which has the gid value
    :param int score_col: index for the column which has the bit score value
    :param int num_lines: number of lines after which a status message is
        logged; defaults to :data:`NUM_LINES`

    :return dict: dictionary for hits
    """

    hits = {}

    for idx, line in enumerate(f_handle):
        if line.startswith('#'):
            continue

        if (idx + 1) % num_lines == 0:
            LOG.info("Parsed %d lines", idx + 1)

        cols = line.strip().split('\t')

        ko_id = cols[0].strip()
        if ko_id in hits:
            continue

        gene_id = cols[gid_col].split('|')[1]

        score = cols[score_col]

        hits[ko_id] = (int(gene_id), float(score))

    return hits


@deprecated
def _parse_gi_taxa_table(f_handle, hits, num_lines=NUM_LINES):
    """
    .. deprecated:: 0.1.13

    Integrates hits dictionary with taxonomic data, parsing the gi taxa table
    from NCBI. Taxon IDs from NCBI are the same as Uniprot taxonomy.

    ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip

    :param file f_handle: the handle for the taxa table from NCBI
    :param dict hits: dictionary with hits data returned by
        :func:`parse_blast_tab`.
    :param int num_lines: number of lines after which a status message is
        logged; defaults to :data:`NUM_LINES`

    :return: dictionary which has ko_idx as key and a namedtuple which has
        ko_id, gi, score, taxon_id.
    """
    gi_ko = {}

    for ko_id, (gene_id, score) in hits.iteritems():
        try:
            gi_ko[gene_id].append(ko_id)
        except KeyError:
            gi_ko[gene_id] = [ko_id]

    # print sum(len(v) for v in gi_ko.itervalues())

    blast_hit = collections.namedtuple('BlastHit', "ko_id gi score taxon_id")

    gi_taxa_dict = {}

    for idx, line in enumerate(f_handle):

        if (idx + 1) % num_lines == 0:
            LOG.info("Parsed %d lines, number of found KO: %d",
                     idx + 1, len(gi_taxa_dict))

        gene_id, taxon_id = line.strip().split()
        gene_id, taxon_id = int(gene_id), int(taxon_id)

        if gene_id in gi_ko:
            for ko_id in gi_ko[gene_id]:
                hit = blast_hit(
                    ko_id=ko_id,
                    gi=gene_id,
                    score=hits[ko_id][1],
                    taxon_id=taxon_id
                )
                # if ko_id in gi_taxa_dict:
                #     print "!"*30
                gi_taxa_dict[ko_id] = hit

    return gi_taxa_dict


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
        ret_col = range(12)

    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading BLAST results from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    for line in file_handle:
        if line.startswith('#'):
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


def parse_uniprot_blast(file_handle, bitscore=40, db='UNIPROT-SP', dbq=10,
                        name_func=None):
    """
    .. versionadded:: 0.1.12

    .. versionchanged:: 0.1.13
        added *name_func* argument

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

    Yields:
        Annotation: instances of :class:`mgkit.io.gff.Annotation` instance of
        each BLAST hit.
    """

    if name_func is None:
        name_func = lambda x: x.split('|')[1]

    # the second function extract the Uniprot ID from the sequence header
    value_funcs = (
        str,
        name_func,
        float,
        int,
        int,
        float
    )

    for seq_id, hit in parse_blast_tab(file_handle, value_funcs=value_funcs):
        if hit[-1] < bitscore:
            continue

        yield gff.from_nuc_blast(hit, db=db, dbq=dbq, to_nuc=False)


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

    for uid, hits in uidmap.iteritems():
        # returns the hit with the max bitscore and max identity
        yield uid, hits


def parse_gi_taxa_table(file_handle, gids=None, num_lines=NUM_LINES):
    """
    .. versionadded:: 0.1.13

    Parses the taxonomy files from the `ncbi ftp
    <ftp://ftp.ncbi.nih.gov/pub/taxonomy/>`_; the file names are
    `gi_taxid_prot*` for the protein db and `gi_taxid_nucl*` for the nucleotide
    db. It contains two columns: the first is the GID and the second is its
    correspondent taxonomy ID.

    Arguments:
        file_handle (str, file): file name or open file handle
        gids (None, list): if it's not `None` only the GIDs included in the
            passed `gids` list will be returned
        num_lines (None, int): number of which a message is logged. If None,
            no message is logged

    Yields:
        tuple: the first element is the GID and the second is the taxonomy ID,
        converted into an integer.

    """
    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')

    if gids is not None:
        gids = set(gids)

    for idx, line in enumerate(file_handle):

        if (num_lines is not None) and ((idx + 1) % num_lines == 0):
            LOG.info("Parsed %d lines", idx + 1)

        gid, taxon_id = line.strip().split('\t')

        if (gids is not None) and (gid not in gids):
            continue

        yield gid, int(taxon_id)
