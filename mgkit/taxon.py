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
from .utils.common import deprecated


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

MISPELLED_TAXA = {
    'bacterioidetes': 'bacteroidetes',
    'bacteriodales': 'bacteroidales'
}
"A few mispelled taxa were present"


UniprotTaxonTuple = collections.namedtuple(
    'UniprotTaxonTuple',
    ('taxon_id', 's_name', 'c_name', 'rank', 'lineage', 'parent_id')
)
"""
A substitute for UniprotTaxon, to be tested
"""


def parse_uniprot_taxon(line):
    """
    .. versionchanged:: 0.1.13
        now accepts empty scientific names, for root taxa

    Parses a Uniprot taxonomy file (tab delimited) line and returns a
    UniprotTaxonTuple instance
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


class UniprotTaxon(object):
    """
    .. deprecated:: 0.1.2

    Class that defines a taxon in Uniprot

    Uses slots to save memory, so no attribute can be attached at instances at
    runtime. To be pickled, __getstate__ and __setstate__ are defined

    At the moment the attributes of a Taxon are:

    * taxon_id (int): id used in Uniprot
    * s_name (string): scientific name
    * c_name (string): common name
    * rank (string): taxon rank (e.g. genus, phylum, etc.)
    * lineage (tuple of strings): the entire lineage that bring to the taxon;
      it's a tuple of strings in the same order as appears in Uniprot
    * parent_id (int): id of the parent taxon

    .. todo::

        get rid of __getstate__ and __setstate__ and test with already pickled
        data

    """

    taxon_id = None
    "Id used in Uniprot (int)"

    s_name = None
    "Scientific name (str)"

    c_name = None
    "Common name (str)"

    rank = None
    "Taxon rank (str)"

    lineage = None
    """
    The entire lineage that bring to the taxon; it's a tuple of strings in the
    same order as appears in Uniprot (tuple of strings)
    """

    parent_id = None
    "Id of the parent taxon (int)"

    @deprecated
    def __init__(self, line=None, **kwd):
        """
        .. deprecated:: 0.1.2

        The method accept any arbitrary keyword arguments, but only accepted
        one will be used

        >>> a="7898\\t9ACTI\\tActinopterygii\\t\\t\\tOsteichthyes; bony fishes; fish; fishes; ray-finned fishes\\t\\tClass\\tEukaryota    ; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi\\t117571"
        >>> tmp1=UniprotTaxon(a)
        >>> tmp1.lineage
        ('eukaryota    ', 'metazoa', 'chordata', 'craniata', 'vertebrata', 'euteleostomi')
        >>> tmp1.taxon_id
        7898
        >>> tmp1.parent_id
        117571
        >>> tmp1.rank
        'class'
        >>> tmp2=UniprotTaxon(**{'lineage': ('eukaryota    ', 'metazoa', 'chordata', 'craniata', 'vertebrata', 'euteleostomi'), 'c_name': '', 'rank': 'class', 'parent_id': 117571, 's_name': 'actinopterygii', 'taxon_id': 7898})
        >>> tmp1.rank == tmp2.rank
        True
        >>> tmp1.lineage == tmp2.lineage
        True
        >>> tmp1.parent_id == tmp2.parent_id
        True
        >>> tmp1.taxon_id == tmp2.taxon_id
        True
        """
        if line is None:
            for name, value in kwd.iteritems():
                setattr(self, name, value)
        else:
            line = line.rstrip().split('\t')
            # if len(line) < 10:
            #     LOG.debug(len(line))
            #     print line
            self.taxon_id = int(line[0])
            self.s_name = line[2].lower()
            try:
                self.c_name = line[3].lower() if line[3] else ''
            except IndexError:
                self.c_name = ''
            try:
                self.rank = line[7].lower()
            except IndexError:
                self.rank = None
            try:
                self.lineage = tuple(x.lower() for x in line[8].split('; '))
            except IndexError:
                self.lineage = (None,)
            try:
                self.parent_id = int(line[9])
            except (ValueError, IndexError):
                self.parent_id = None

    def __getstate__(self):
        """
        Used by pickle.dump

        >>> tmp1=UniprotTaxon(**{'lineage': ('eukaryota    ', 'metazoa', 'chordata', 'craniata', 'vertebrata', 'euteleostomi'), 'c_name': '', 'rank': 'class', 'parent_id': 117571, 's_name': 'actinopterygii', 'taxon_id': 7898})
        >>> tmp1.__getstate__()
        {'lineage': ('eukaryota    ', 'metazoa', 'chordata', 'craniata', 'vertebrata', 'euteleostomi'), 'c_name': '', 'rank': 'class', 'parent_id': 117571, 's_name': 'actinopterygii', 'taxon_id': 7898}
        """
        return dict((x, getattr(self, x)) for x in self.__dict__)

    def __setstate__(self, state):
        """
        Used by pickle.load
        """
        for name, value, in state.iteritems():
            setattr(self, name, value)

    def __repr__(self):
        return str(self)

    def __str__(self):
        """
        String representation for the class instance

        >>> UniprotTaxon(**{'lineage': ('eukaryota    ', 'metazoa', 'chordata', 'craniata', 'vertebrata', 'euteleostomi'), 'c_name': '', 'rank': 'class', 'parent_id': 117571, 's_name': 'actinopterygii', 'taxon_id': 7898})
        actinopterygii (7898) - class
        """
        return "{0} ({1}) - {2}".format(self.s_name, self.taxon_id, self.rank)


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

    def read_taxonomy(self, f_handle):
        """
        Reads taxonomy from a file handle.
        The file needs to be a tab separated format return by a query on
        Uniprot.

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
                taxon = parse_uniprot_taxon(line)
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
        elif isinstance(taxon_id, UniprotTaxon):
            ranked = taxon_id
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
        elif isinstance(taxon_id, UniprotTaxon):
            taxon = taxon_id
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
        if isinstance(taxon, UniprotTaxon):
            return taxon in self._taxa.itervalues()
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


def reorder_iterable_by_root_taxon(taxa_list, tmap, aggr=None):
    """
    Reorders iterable based on a root taxa map.

    :param iterable taxa_list: any iterable with taxa names
    :param dict tmap: dictionary returned by :func:`group_by_root` with
        only_names=True
    :param func aggr: callable to aggregate results. (e.g. list, pandas.Index,
        ecc.)

    :return: a generator with the ordered taxa or the data that aggr returns

    .. todo::

        * option to output each taxa list for a root taxon sorted
    """

    ord_index = {}

    for taxon in taxa_list:
        try:
            ord_index[get_taxon_root(taxon, tmap)].append(taxon)
        except KeyError:
            ord_index[get_taxon_root(taxon, tmap)] = [taxon]

    if aggr is None:
        return itertools.chain(*ord_index.values())
    else:
        return aggr(itertools.chain(*ord_index.values()))


def get_taxon_root(taxon, tmap, replace_space=False):
    """
    Returns a taxon root, given a taxon name and a taxa map

    :param str taxon: taxon name
    :param tmap: dictionary returned by :func:`group_by_root` with
        only_names=True
    :param bool replace_space: if True replaces '#' with a space before the
        check

    :return: root taxon name

    .. todo::

        make sure that the mispelled taxa are fixed (in hmmr2gff) and that the
        other fix (subclass rank included in the pickled data) is done
    """
    #mispelled, change gff file?
    if taxon in ['bacterioidetes', 'bacteriodales']:
        return 'bacteria'

    #makes sure taxon names with spaces can be used
    if replace_space:
        taxon = taxon.replace('#', ' ')

    for root, taxa in tmap.iteritems():
        if (taxon == root) or (taxon in taxa):
            return root
    raise ValueError("Taxon '{0}' not in mapping".format(taxon))


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
