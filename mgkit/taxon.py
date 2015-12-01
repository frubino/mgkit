"""
This module gives access to Uniprot taxonomy data. It also defines classes
to filter, order and group data by taxa
"""

import logging
import gzip
import cPickle
import itertools
import collections
import pandas
from .io import open_file
# from .utils.common import deprecated


LOG = logging.getLogger(__name__)

ALGAE = {
    'haptophyceae': 2830,
    'chlorophyta': 3041,
    'stramenopiles': 33634,
    'cryptophyta': 3027,
    'rhodophyta': 2763,
}

BACTERIA = 2
ARCHAEA = 2157
FUNGI = 4751

PROTISTS = {
    #'apicomplexa', alveolata
    #'ciliophora', alveolata
    'alveolata': 33630,
    'amoebozoa': 554915,
    'apusozoa': 554296,
    'breviatea': 1401294,
    'centroheliozoa': 193537,
    'choanoflagellida': 28009,
    'diplomonadida': 5738,
    'euglenozoa': 33682,
    'formicata': 207245,
    'heterolobosea': 5752,
    'jakobida': 556282,
    'malawimonadidae': 136087,
    'oxymonadida': 66288,
    'parabasalia': 5719,
    'rhizaria': 543769,
}

PLANTS = {
    'streptophyta': 35493,
}

TAXON_ROOTS = (
    'archaea',  # 2157
    'bacteria',  # 2
    'fungi',  # 4751
    'metazoa',
    'environmental samples',
    'viruses',
    'viroids'
    'eukaryota',
    'other sequences',
    'unidentified'
) + tuple(PROTISTS) + tuple(PLANTS) + tuple(ALGAE)
"Root taxa used in analysis and filtering"

TAXON_RANKS = (
    'superkingdom',
    'kingdom',
    'phylum',
    'class',
    'subclass',
    'order',
    'family',
    'genus',
    'species'
)
"Taxonomy ranks included in the pickled data"

UniprotTaxonTuple = collections.namedtuple(
    'UniprotTaxonTuple',
    ('taxon_id', 's_name', 'c_name', 'rank', 'lineage', 'parent_id')
)
"""
A representation of a Uniprot Taxon
"""


def parse_uniprot_taxon(line, light=True):
    """
    .. versionchanged:: 0.1.13
        now accepts empty scientific names, for root taxa

    .. versionchanged:: 0.2.1
        added *light* parameter

    Parses a Uniprot taxonomy file (tab delimited) line and returns a
    UniprotTaxonTuple instance. If *light* is True, lineage is not stored to
    decrease the memory usage. This is now the default.
    """
    line = line.rstrip().split('\t')
    taxon_id = int(line[0])
    try:
        s_name = line[2].lower()
    except IndexError:
        s_name = ''
    try:
        c_name = line[3].lower() if line[3] else ''
    except IndexError:
        c_name = ''
    try:
        rank = line[7].lower()
    except IndexError:
        rank = None
    try:
        if light:
            raise IndexError
        lineage = tuple(x.lower() for x in line[8].split('; '))
    except IndexError:
        lineage = (None,)
    try:
        parent_id = int(line[9])
    except (ValueError, IndexError):
        parent_id = None
    return UniprotTaxonTuple(
        taxon_id,
        s_name,
        c_name,
        rank,
        lineage,
        parent_id
    )


class UniprotTaxonomy(object):
    """
    Class that contains the whole Uniprot taxonomy. Defines some methods to
    easy access of taxonomy.

    Defines:

    * methods to load taxonomy from a pickle file or a generic file handle
    * can be iterated over and returns a generator its UniprotTaxon instances
    * can be used as a dictionary, in which the key is a taxon_id and the value
      is its UniprotTaxon instance
    """
    def __init__(self, fname=None):
        """
        :param str fname: file name of the pickled data.
        """
        self._taxa = {}
        self._name_map = {}
        if fname:
            self.load_data(fname)

    def read_taxonomy(self, f_handle, light=True):
        """
        .. versionchanged:: 0.2.1
            added *light* parameter

        Reads taxonomy from a file handle.
        The file needs to be a tab separated format return by a query on
        Uniprot.  If *light* is True, lineage is not stored to decrease the
        memory usage. This is now the default.

        New taxa will be added, duplicated taxa will be skipped.

        :param handle f_handle: file handle of the taxonomy file.
        """
        try:
            LOG.info("Loading taxonomy from file %s", f_handle.name)
        except AttributeError:
            pass

        f_handle.readline()

        for line in f_handle:
            try:
                # taxon = UniprotTaxon(line=line)
                taxon = parse_uniprot_taxon(line, light=light)
            except Exception:
                # print line.strip().split('\t')
                LOG.debug("skipped line: %s", line)
                continue
            if taxon.taxon_id in self:
                continue
            self[taxon.taxon_id] = taxon

    def gen_name_map(self):
        """
        Generate a name map, where to each scientific name in the taxonomy an
        id is associated.
        """
        LOG.debug('Generate name map')
        self._name_map = {}
        for taxon_obj in self:
            try:
                self._name_map[taxon_obj.s_name].append(taxon_obj.taxon_id)
            except KeyError:
                self._name_map[taxon_obj.s_name] = [taxon_obj.taxon_id]

    def load_data(self, file_handle):
        """
        .. versionchanged:: 0.1.13
            now accepts file handles and compressed files (if file names)

        Loads pickled data from file name "file_handle", accept gzipped data

        :param str file_handle: file name (or file handle) of the pickled data
        """

        if isinstance(file_handle, str):
            file_handle = open_file(file_handle)

        LOG.info("Loading taxonomy from file %s", file_handle.name)

        self._taxa = cPickle.load(file_handle)

    def save_data(self, fname):
        """
        Saves taxonomy data to file name "fname", can write gzipped data if
        fname ends with ".gz"

        :param str fname: file name to which save the instance data
        """
        LOG.info("Saving taxonomy to file %s", fname)
        if fname.endswith('.gz'):
            f_handle = gzip.open(fname, 'w', compresslevel=4)
        else:
            f_handle = open(fname, 'w')
        cPickle.dump(self._taxa, f_handle, -1)

    def find_by_name(self, s_name):
        """
        Returns the taxon IDs associated with the scientific name provided

        :param str s_name: the scientific name

        :return list: a reference to the list of IDs that have for `s_name`
        """
        if not self._name_map:
            self.gen_name_map()
        return self._name_map[s_name]

    def is_ancestor(self, leaf_id, anc_ids):
        """
        .. versionchanged:: 0.1.13
            now uses :func:`is_ancestor` and changed behavior

        Checks if a taxon is the leaf of another one, or a list of taxa.

        :param int leaf_id: leaf taxon id
        :param int anc_ids: ancestor taxon id(s)

        :return bool: True if the ancestor taxon is in the leaf taxon lineage

        """

        if isinstance(anc_ids, int):
            anc_ids = [anc_ids]

        for anc_id in anc_ids:
            if is_ancestor(self, leaf_id, anc_id):
                return True
        return False

    def get_ranked_taxon(self, taxon_id, rank=None, ranks=TAXON_RANKS, roots=False):
        """
        .. versionchanged:: 0.1.13
            added *roots* argument

        Traverse the branch of which the *taxon* argument is the leaf backward,
        to get the specific rank to which the *taxon* belongs to.

        .. warning::

            the *roots* options is kept for backward compatibility and should be
            be set to *False*

        :param taxon_id: id of the taxon or instance of :class:`UniprotTaxon`
        :param str rank: string that specify the rank, if None, the first valid
            rank will be searched. (i.e. the first with a value different from '')
        :param ranks: tuple of all taxonomy ranks, default to the default module
            value
        :param bool roots: if True, uses :data:`TAXON_ROOTS` to solve the root
            taxa
        :return: instance of :class:`UniprotTaxon` for the rank found.
        """

        if isinstance(taxon_id, int):
            ranked = self[taxon_id]
        else:
            raise ValueError("Not a valid taxon_id: {0}".format(taxon_id))

        while True:
            if (rank is None) and (ranked.rank != ''):
                break
            elif (ranked.rank == rank):
                break
            elif (ranked.rank in ranks) and \
                    (ranks.index(ranked.rank) < ranks.index(rank)):
                break
            elif ranked.parent_id is None:
                break
            #kept only for backward compatibility.
            #needs a check to other scripts
            elif roots and (ranked.s_name in TAXON_ROOTS):
                break

            ranked = self[ranked.parent_id]

        return ranked

    def get_name_map(self):
        "Returns a taxon_id->s_name dictionary"

        return dict(
            (taxon.taxon_id, taxon.s_name)
            for taxon in self
        )

    def get_taxon_root(self, taxon_id, roots=TAXON_ROOTS):
        """
        Given a :class:`UniprotTaxon` instance and the associated taxonomy, returns
        the correct root taxon the supplied taxon belongs to.

        :param int taxon_id: taxon id or instance of :class:`UniprotTaxon`
        :param iterable roots: list of root taxa scientific names

        :return: root :class:`UniprotTaxon` instance
        """

        if isinstance(taxon_id, int):
            taxon = self[taxon_id]
        else:
            raise ValueError("Not a valid taxon_id: {0}".format(taxon_id))

        while taxon.s_name not in roots:
            taxon = self[taxon.parent_id]

        return taxon

    def __getitem__(self, taxon_id):
        """
        Defines dictionary behavior. Key is a taxon_id, the returned value is a
        UniprotTaxon instance
        """
        return self._taxa[taxon_id]

    def __setitem__(self, taxon_id, taxon_obj):
        self._taxa[taxon_id] = taxon_obj

    def __delitem__(self, taxon_id):
        del self._taxa[taxon_id]

    def __iter__(self):
        """
        Defines iterable behavior. Returns a generator for UniprotTaxon instances
        """
        for taxon in self._taxa.itervalues():
            yield taxon

    def __len__(self):
        """
        Returns the number of taxa contained
        """
        return len(self._taxa)

    def __contains__(self, taxon):
        """
        Returns True if the taxon is in the taxonomy

        Accepts an int (check for taxon_id) or an instance of UniprotTaxon
        """
        if isinstance(taxon, int):
            return taxon in self._taxa
        return False


def group_by_root(taxa, roots=TAXON_ROOTS, only_names=False, replace_space='#'):
    """
    Returns a dictionary containing as keys the root taxa and as values a list
    of taxa belonging to that group (checks lineage attribute of
    :class:`UniprotTaxon` instances)

    :param iterable taxa: iterable of :class:`UniprotTaxon` instances
    :param roots: root taxa to be used to construct the dictionary, defaults to
        the 5 defined in :data:`TAXON_ROOTS`
    :param bool only_names: boolean that specify what data type to return for
        values, if True returns strings (taxa names), if False
        :class:`UniprotTaxon` instances
    :param str replace_space: if only_names is True replaces spaces with the
        choosen character

    :return: a dictionary root->[taxa_ids]
    """
    groups = dict(
        (root, []) for root in roots
    )
    if not only_names:
        taxa.gen_name_map()
        for root in roots:
            groups[root].append(taxa.find_by_name(root)[0])

    for taxon in taxa:
        for root in roots:
            if root in taxon.lineage:
                groups[root].append(
                    taxon.s_name.replace(
                        ' ', replace_space) if only_names else taxon.taxon_id
                )
                break

    return groups


def load_taxonomy_map(taxon_file):
    "Loads taxonomy from file and return a map root_taxon->taxa"
    taxonomy = UniprotTaxonomy(taxon_file)
    return group_by_root(taxonomy, only_names=True)


def get_ancestor_map(leaf_ids, anc_ids, taxonomy):
    """
    This function returns a dictionary where every leaf taxon is associated
    with the right ancestors in anc_ids

    ex. {clostridium: [bacteria, clostridia]}
    """
    anc_map = {}

    for leaf_id in leaf_ids:
        for anc_id in anc_ids:
            if taxonomy.is_ancestor(leaf_id, anc_id) or (leaf_id == anc_id):
                try:
                    anc_map[leaf_id].add(anc_id)
                except KeyError:
                    anc_map[leaf_id] = set()
                    anc_map[leaf_id].add(anc_id)

    return anc_map


def is_ancestor(taxonomy, taxon_id, anc_id):
    """
    .. versionchanged:: 0.1.16
        if a *taxon_id* raises a KeyError, False is returned

    Determine if the given taxon id (taxon_id) has anc_id as ancestor.

    :param :class:`UniprotTaxonomy` taxonomy: taxonomy used to test
    :param int taxon_id: leaf taxon to test
    :param int anc_id: ancestor taxon to test against

    :return bool: True if anc_id is an ancestor of taxon_id or their the same
    """
    if taxon_id == anc_id:
        return True

    while True:
        # if a taxon_id is not found, ancestry is invalid anyway, so it returns
        # False
        try:
            taxon_id = taxonomy[taxon_id].parent_id
        except KeyError:
            return False
        # print taxon_id, anc_id
        if anc_id == taxon_id:
            return True
        if taxon_id is None:
            return False


class NoLcaFound(Exception):
    """
    .. versionadded:: 0.1.13

    Raised if no lowest common ancestor can be found in the taxonomy
    """
    pass


def last_common_ancestor(taxonomy, taxon_id1, taxon_id2):
    """
    .. versionadded:: 0.1.13

    Finds the last common ancestor of two taxon IDs. An alias to this function
    is in the same module, called *lowest_common_ancestor* for compatibility.

    Arguments:
        taxonomy: :class:`UniprotTaxonomy` instance used to test
        taxon_id1 (int): first taxon ID
        taxon_id2 (int): second taxon ID

    Raturns:
        int: taxon ID of the lowest common ancestor

    Raises:
        NoLcaFound: if no common ancestor can be found
    """
    lca_id = taxon_id1

    while True:
        if is_ancestor(taxonomy, taxon_id2, lca_id):
            break
        else:
            lca_id = taxonomy[lca_id].parent_id
            if lca_id is None:
                raise NoLcaFound('No common ancestry')

    return lca_id

lowest_common_ancestor = last_common_ancestor


def distance_taxa_ancestor(taxonomy, taxon_id, anc_id):
    """
    .. versionadded:: 0.1.16

    Function to calculate the distance between a taxon and the given ancestor

    The distance is equal to the number of step in the taxonomy taken to arrive
    at the ancestor.

    Arguments:
        taxonomy: :class:`UniprotTaxonomy` instance
        taxon_id (int): taxonomic identifier
        anc_id (int): taxonomic identifier of the ancestor

    Raturns:
        int: distance between *taxon_id* and it ancestor *anc_id*

    """
    distance = 0

    while True:
        if anc_id == taxon_id:
            break
        taxon_id = taxonomy[taxon_id].parent_id
        distance += 1

    return distance


def distance_two_taxa(taxonomy, taxon_id1, taxon_id2):
    """
    .. versionadded:: 0.1.16

    Calculate the distance between two taxa. The distance is equal to the sum
    steps it takes to traverse the taxonomy until their last common ancestor.

    Arguments:
        taxonomy: :class:`UniprotTaxonomy` instance
        taxon_id1 (int): taxonomic identifier of first taxon
        taxon_id2 (int): taxonomic identifier of second taxon

    Raturns:
        int: distance between *taxon_id1* and *taxon_id2*
    """
    lca_id = last_common_ancestor(taxonomy, taxon_id1, taxon_id2)

    return sum(
        distance_taxa_ancestor(taxonomy, taxon_id, lca_id)
        for taxon_id in (taxon_id1, taxon_id2)
    )


def taxa_distance_matrix(taxonomy, taxon_ids):
    """
    .. versionadded:: 0.1.16

    Given a list of taxonomic identifiers, returns a distance matrix in a
    pairwise manner by using :func:`distance_two_taxa` on all possible
    two element combinations of *taxon_ids*.

    Arguments:
        taxonomy: :class:`UniprotTaxonomy` instance
        taxon_ids (iterable): list taxonomic identifiers

    Returns:
        pandas.DataFrame: matrix with the pairwise distances of all *taxon_ids*
    """
    matrix = pandas.DataFrame(index=taxon_ids, columns=taxon_ids).fillna(0)

    for taxon_id1, taxon_id2 in itertools.combinations(taxon_ids, 2):
        distance = distance_two_taxa(taxonomy, taxon_id1, taxon_id2)
        matrix.loc[taxon_id1, taxon_id2] = distance
        matrix.loc[taxon_id2, taxon_id1] = distance

    return matrix


def get_lineage(taxonomy, taxon_id, names=False):
    """
    .. versionadded:: 0.2.1

    Returns the lineage of a taxon_id, as a list of taxon_id or taxa names

    Arguments:
        taxonomy: a :class:`UniprotTaxonomy` instance
        taxon_id (int): taxon_id whose lineage to return
        names (bool): if True, the returned list contains the names of the taxa
            instead of the taxon_id

    Returns:
        list: lineage of the taxon_id, the elements are `int` if names is False,
        and `str` when names is True
    """
    lineage = []

    while True:
        taxon_id = taxonomy[taxon_id].parent_id
        if taxon_id is None:
            break
        lineage.append(taxon_id)

    if names:
        lineage = [taxonomy[x].s_name for x in lineage]

    return list(reversed(lineage))
