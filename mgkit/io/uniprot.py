"""
.. versionadded:: 0.1.13

Uniprot file formats
"""
import sys
import logging
from . import open_file

LOG = logging.getLogger(__name__)

NUM_LINES = 10 ** 7


MAPPINGS = {
    'taxonomy': 'NCBI_TaxID',
    'eggnog': 'eggNOG',
    'ko': 'KO',
    'kegg': 'KEGG',
    'biocyc': 'BioCyc',
    'unipathway': 'UniPathway',
    'embl': 'EMBL',
    'embl_cds': 'EMBL-CDS',
    'gi': 'GI',
    'string': 'STRING'
}
"""
Some of the mappings contained in the idmapping.dat.gz
"""


def parse_uniprot_mappings(file_handle, gene_ids=None, mappings=None,
                           num_lines=NUM_LINES):
    """
    Parses a Uniprot mapping `file <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_,
    returning a generator with the mappings.

    Arguments:
        file_handle (str, file): file name or open file handle
        gene_ids (None, set): if not None, the returned mappings are for the
            gene IDs specified
        mappings (None, set): mappings to be returned
        num_lines (None, int): number of which a message is logged. If None,
            no message is logged

    Yields:
        tuple: the first element is the gene ID, the second is the mapping type
        and third element is the mapped ID
    """
    if (sys.version_info[0] == 2) and isinstance(file_handle, unicode):
        file_handle = open_file(file_handle, 'rb')
    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'rb')

    LOG.info(
        "Loading Uniprot Mappings from file (%s)",
        getattr(file_handle, 'name', repr(file_handle))
    )

    if gene_ids is not None:
        gene_ids = set(gene_ids)
        LOG.info("Mappings for %d gene_ids will be returned", len(gene_ids))

    if mappings is not None:
        mappings = set(mappings)
        LOG.info(
            "Mappings to '%s' will be returned",
            ', '.join(mappings)
        )

    for idx, line in enumerate(file_handle):
        line = line.decode('ascii')

        if (num_lines is not None) and ((idx + 1) % num_lines == 0):
            LOG.info("Parsed %d lines", idx + 1)

        gene_id, mapping, map_id = line.strip().split('\t')

        if (gene_ids is not None) and (gene_id not in gene_ids):
            continue

        if (mappings is not None) and (mapping not in mappings):
            continue

        yield gene_id, mapping, map_id

    LOG.info("Read %d lines", idx + 1)


def uniprot_mappings_to_dict(file_handle, gene_ids, mappings, num_lines=None):
    """
    .. versionchanged:: 0.3.4
        added *num_lines*

    Parses a Uniprot mapping `file <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_,
    returning a generator of dictionaries with the mappings requested.

    Arguments:
        file_handle (str, file): file name or open file handle
        gene_ids (None, set): if not None, the returned mappings are for the
            gene IDs specified
        mappings (None, set): mappings to be returned
        num_lines (int, None): passed to :func:`parse_uniprot_mappings`

    Yields:
        tuple: the first element is the gene ID, the second is a dictionary
        with all the mappings found, the key is the mapping type and the value
        is a list of all mapped IDs
    """
    iterator = parse_uniprot_mappings(
        file_handle,
        gene_ids=gene_ids,
        mappings=mappings,
        num_lines=num_lines
    )

    curr_gene = ''
    curr_maps = {}

    for gene_id, mapping, map_id in iterator:
        if curr_gene == gene_id:
            try:
                curr_maps[mapping].append(map_id)
            except KeyError:
                curr_maps[mapping] = [map_id]
        else:
            if curr_gene != '':
                yield curr_gene, curr_maps
                curr_maps = {}
            curr_gene = gene_id
            curr_maps[mapping] = [map_id]
    else:
        yield curr_gene, curr_maps
