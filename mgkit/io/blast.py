"""
Blast routines and parsers

"""

import logging
import collections

NUM_LINES = 10**6

LOG = logging.getLogger(__name__)


def parse_blast_tab(f_handle, gid_col=1, score_col=11, num_lines=NUM_LINES):
    """
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

        gi = cols[gid_col].split('|')[1]

        score = cols[score_col]

        hits[ko_id] = (int(gi), float(score))

    return hits


def parse_gi_taxa_table(f_handle, hits, num_lines=NUM_LINES):
    """
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

    for ko_id, (gi, score) in hits.iteritems():
        try:
            gi_ko[gi].append(ko_id)
        except KeyError:
            gi_ko[gi] = [ko_id]

    # print sum(len(v) for v in gi_ko.itervalues())

    blast_hit = collections.namedtuple('BlastHit', "ko_id gi score taxon_id")

    gi_taxa_dict = {}

    for idx, line in enumerate(f_handle):

        if (idx + 1) % num_lines == 0:
            LOG.info("Parsed %d lines, number of found KO: %d",
                     idx + 1, len(gi_taxa_dict))

        gi, taxon_id = line.strip().split()
        gi, taxon_id = int(gi), int(taxon_id)

        if gi in gi_ko:
            for ko_id in gi_ko[gi]:
                hit = blast_hit(
                    ko_id=ko_id,
                    gi=gi,
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
        #skips if the ID is not in the Uniprot Taxonomy
        annotation.attributes.blast_taxon = taxonomy[hit.taxon_id].s_name
        annotation.attributes.blast_taxon_idx = hit.taxon_id
    except KeyError:
        return
