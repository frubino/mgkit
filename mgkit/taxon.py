"""
This module gives access to Uniprot taxonomy data. It also defines classes
to filter, order and group data by taxa
"""
from builtins import object
from functools import reduce
import logging
import sys
if sys.version_info[0] == 2:
    import cPickle as pickle
else:
    import pickle
import itertools
import collections
from future.utils import viewitems, viewvalues
from .io import open_file
from . import DependencyError


LOG = logging.getLogger(__name__)

ALGAE = {
    'haptophyceae': 2830,
    'chlorophyta': 3041,
    'stramenopiles': 33634,
    'cryptophyta': 3027,
    'rhodophyta': 2763,
}

# superkingdoms
BACTERIA = 2
ARCHAEA = 2157
VIRUS = 10239

CELLULAR_ORGANISMS = 131567

# kingdoms
FUNGI = 4751
VIRIDIPLANTAE = 33090
METAZOA = 33208
EUKARYOTA = 2759

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

TaxonTuple = collections.namedtuple(
    'UniprotTaxonTuple',
    ('taxon_id', 's_name', 'c_name', 'rank', 'lineage', 'parent_id')
)
"""
A representation of a Uniprot Taxon
"""

# To be deprecated to avoid confusion
UniprotTaxonTuple = TaxonTuple


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
    line = line.decode('ascii')
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


def parse_ncbi_taxonomy_merged_file(file_handle):
    """
    .. versionadded:: 0.2.3

    Parses the *merged.dmp* file where the merged taxon_id are stored. Available
    at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

    Arguments:
        file_handle (str, file): file name or handle to the file

    Returns:
        dict: dictionary with merged_id -> taxon_id
    """
    file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading NCBI taxonomy merged file %s",
        getattr(file_handle, 'name', repr(file_handle))
    )

    merged_taxa = {}

    for line in file_handle:
        line = line.decode('ascii')
        merged_id, taxon_id = [col for col in line.strip().split('\t') if col != '|']
        merged_taxa[int(merged_id)] = int(taxon_id)

    return merged_taxa


def parse_ncbi_taxonomy_names_file(file_handle, name_classes=('scientific name', 'common name')):
    """
    .. versionadded:: 0.2.3

    Parses the *names.dmp* file where the names associated to a taxon_id are
    stored. Available at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

    Arguments:
        file_handle (str, file): file name or handle to the file
        name_classes (tuple): name classes to save, only the scientific and
            common name are stored

    Returns:
        dict: dictionary with merged_id -> taxon_id
    """

    file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading NCBI taxonomy names file %s",
        getattr(file_handle, 'name', repr(file_handle))
    )

    taxa_names = {}

    for line in file_handle:
        line = line.decode('ascii')
        taxon_id, taxon_name, uniq_name, name_class = [
            col.strip()
            for col in line.strip().split('\t')
            if col != '|'
        ]
        taxon_id = int(taxon_id)

        if name_class not in name_classes:
            continue

        if taxon_id not in taxa_names:
            taxa_names[taxon_id] = {}

        taxa_names[taxon_id][name_class] = taxon_name

    return taxa_names


def parse_ncbi_taxonomy_nodes_file(file_handle, taxa_names=None):
    """
    .. versionadded:: 0.2.3

    Parses the *nodes.dmp* file where the nodes of the taxonomy are stored.
    Available at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/.

    Arguments:
        file_handle (str, file): file name or handle to the file
        taxa_names (dict): dictionary with the taxa names (returned from
            :func:`parse_ncbi_taxonomy_names_file`)

    Yields:
        TaxonTuple: TaxonTuple instance
    """

    file_handle = open_file(file_handle, 'r')

    LOG.info(
        "Reading NCBI taxonomy nodes file %s",
        getattr(file_handle, 'name', repr(file_handle))
    )

    for line in file_handle:
        line = line.decode('ascii')
        line = [col for col in line.strip().split('\t') if col != '|']
        taxon_id = int(line[0])
        parent_id = int(line[1])
        if parent_id == 1:
            # NCBI uses 1 as the root, but the rest of the functions expect None
            parent_id = None
        rank = line[2].lower()

        s_name = ''
        c_name = ''

        if taxa_names is not None:
            try:
                names = taxa_names[taxon_id]
            except KeyError:
                names = {}
            s_name = names.get('scientific name', '')
            c_name = names.get('common name', '')

        yield TaxonTuple(
            taxon_id,
            s_name,
            c_name,
            rank,
            (None,),  # lineage is not found in the NCBI taxonomy dump
            parent_id
        )


class Taxonomy(object):
    """
    Class that contains the whole Uniprot taxonomy. Defines some methods to
    easy access of taxonomy. Follows the conventions of NCBI Taxonomy.

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
    @staticmethod
    def parse_gtdb_lineage(lineage, sep=';'):
        """
        .. versionadded:: 0.3.3

        Parse a GTDB lineage, one that defines the rank as a single letter,
        followed by **__** for each taxon name. Taxa are separated by semicolon
        by default. Also the **domain** rank is renamed into **superkingdom**
        to allow mixing of taxonomies.

        Returns:
            dict: dictionary with the parsed lineage, which can be passed to
            :meth:`Taxonomy.add_lineage`
        """
        lineage_dict = {}
        ranks = dict(d='superkingdom', f='family', p='phylum', o='order', g='genus', s='species', c='class')
        for taxon_name in lineage.split(sep):
            rank, taxon_name = taxon_name.split('__')
            if not taxon_name:
             continue
            lineage_dict[ranks[rank]] = taxon_name
        return lineage_dict

    def read_from_gtdb_taxonomy(self, file_handle, use_gtdb_name=True, sep='\t'):
        """
        .. versionadded:: 0.3.0

        .. versionchanged:: 0.3.1
            replaced *domain* with *superkingdom* to support *get_lineage*

        Reads a GTDB taxonomy file (tab separated genome_id/taxonomy) and
        populate the taxonomy instance. The method also return a dictionary of
        genome_id -> taxon_id.

        Arguments:
            file_handle (file): file with the taxonomy
            use_gtdb_name (bool): if True, the names are kept as-is in the
                *s_name* attribute of :class:`TaxonTuple` and the
                "cleaned" version in *c_name* (e.g. f__Ammonifexaceae ->
                Ammonifexaceae). If False, the values are switched
            sep (str): separator between the columns of the file

        Returns:
            dict: dictionary of genome_id -> taxon_id, reflecting the created
            taxonomy

        .. note::

            the taxon_id are generated, so there's no guarantee they will be
            the same in a successive execution

        """
        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'r')

        LOG.info(
            "Reading GTDB taxonomy from file (%s)",
            getattr(file_handle, 'name', repr(file_handle))
        )

        genome_ids = {}

        # ranks used in GTDB
        ranks = {
            'g': 'genus',
            'f': 'family',
            'p': 'phylum',
            'o': 'order',
            'd': 'superkingdom',
            's': 'species',
            'c': 'class'
        }

        taxon_ids = {}
        count = 1
        for line in file_handle:
            line = line.decode('ascii')
            # expecting the table to be the exported file from GTDB
            try:
                # print line.strip().split(sep)
                genome_id, line = line.strip().split(sep)
                # in case there's no gtdb taxonomy
            except ValueError:
                continue

            line = [
                taxon_name.strip()
                for taxon_name in line.split(';')
                if len(taxon_name) > 3
            ]

            if not line:
                continue

            # if no parent, use None as the NCBI taxonomy
            parent_id = None
            for taxon_name in line:

                # register a taxo_id for unknown taxon
                if taxon_name not in taxon_ids:
                    taxon_ids[taxon_name] = count
                    count += 1

                taxon_id = taxon_ids[taxon_name]
                try:
                    rank = ranks[taxon_name[0]]
                except KeyError:
                    raise KeyError

                # keep the full gtdb name in the common name
                # cut the scientific name to remove the rank information
                common_name = taxon_name
                taxon_name = taxon_name[3:]
                # but if th use_gtdb_name option is True, switch them
                if use_gtdb_name:
                    taxon_name, common_name = common_name, taxon_name

                self[taxon_id] = TaxonTuple(
                    taxon_id,
                    taxon_name,
                    common_name,
                    rank,
                    (None,),
                    parent_id
                )
                parent_id = taxon_id
            genome_ids[genome_id] = taxon_id

        return genome_ids

    def read_from_ncbi_dump(self, nodes_file, names_file=None, merged_file=None):
        """
            .. versionadded:: 0.2.3

            Uses the *nodes.dmp* and optionally *names.dmp*, *merged.dmp* files
            from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ to populate the
            taxonomy.

            Arguments:
                nodes_file (str, file): file name or handle to the file
                names_file (str, file, None): file name or handle to the file,
                    if None, names won't be added to the taxa
                merged_file (str, file, None): file name or handle to the file,
                    if None, pointers to merged taxa won't be added
            """

        LOG.info("Reading NCBI taxonomy from file")

        if names_file is not None:
            taxa_names = parse_ncbi_taxonomy_names_file(names_file)
        else:
            taxa_names = None

        for taxon in parse_ncbi_taxonomy_nodes_file(nodes_file, taxa_names):
            self[taxon.taxon_id] = taxon

        if merged_file is not None:
            merged_taxa = parse_ncbi_taxonomy_merged_file(merged_file)

            for merged_id, taxon_id in viewitems(merged_taxa):
                self[merged_id] = self[taxon_id]

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
        .. versionchanged:: 0.2.3
            names are stored in the mapping as lowercase

        Generate a name map, where to each scientific name in the taxonomy an
        id is associated.
        """
        LOG.debug('Generate name map')
        self._name_map = {}
        for taxon_obj in self:
            s_name = taxon_obj.s_name.lower()
            try:
                self._name_map[s_name].append(taxon_obj.taxon_id)
            except KeyError:
                self._name_map[s_name] = [taxon_obj.taxon_id]

    def load_data(self, file_handle):
        """
        .. versionchanged:: 0.2.3
            now can use read *msgpack* serialised files

        .. versionchanged:: 0.1.13
            now accepts file handles and compressed files (if file names)

        Loads serialised data from file name "file_handle" and accepts
        compressed files.

        if the *.msgpack* string is found in the file name, the *msgpack*
        package is used instead of pickle

        Arguments:
            file_handle (str, file): file name to which save the instance data

        Raises:
            DependencyError: if the file name contains *.msgpack* and the
            package is not installed
        """

        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'rb')

        LOG.info("Loading taxonomy from file %s", file_handle.name)

        if '.msgpack' in file_handle.name:
            try:
                import msgpack
            except ImportError:
                raise DependencyError('msgpack')
            merged = []
            for taxon_id, taxon in msgpack.Unpacker(file_handle, use_list=False, raw=False):
                taxon = TaxonTuple(*taxon)
                # if it's a merged taxon_id keep it and don't add it
                if taxon_id != taxon.taxon_id:
                    merged.append((taxon_id, taxon.taxon_id))
                else:
                    self[taxon_id] = taxon
            # the list of merged ids is iterated over and pointers are added to
            # the correct ids
            for merged_id, taxon_id in merged:
                self[merged_id] = self[taxon_id]

        else:
            self._taxa = pickle.load(file_handle)

    def save_data(self, file_handle):
        """
        .. versionchanged:: 0.2.3
            now can use *msgpack* to serialise

        Saves taxonomy data to a file handle or file name, can write compressed
        data if the file ends with ".gz", ".bz2"

        if the *.msgpack* string is found in the file name, the *msgpack*
        package is used instead of pickle

        Arguments:
            file_handle (str, file): file name to which save the instance data

        Raises:
            DependencyError: if the file name contains *.msgpack* and the
            package is not installed
        """

        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'wb')

        LOG.info("Saving taxonomy to file %s", file_handle.name)

        if '.msgpack' in file_handle.name:
            try:
                import msgpack
            except ImportError:
                raise DependencyError('msgpack')
            for taxon_id, taxon in viewitems(self._taxa):
                file_handle.write(
                    msgpack.packb((taxon_id, taxon))
                )
        else:
            pickle.dump(self._taxa, file_handle, -1)

    def find_by_name(self, s_name, rank=None, strict=True):
        """
        .. versionchanged:: 0.2.3
            the search is now case insensitive

        .. versionchanged:: 0.3.1
            added *rank* and *strict* parameter

        Returns the taxon IDs associated with the scientific name provided

        Arguments:
            s_name (str): the scientific name
            rank (str, None): return only a taxon_id of a specific rank
            strict (book): if True and *rank* is not None, KeyError will be
                raised if multiple taxa have the same name and rank

        Returns:
            list: a reference to the list of IDs that have for *s_name*, if
            *rank* is None. If *rank* is not None and one taxon is found, its
            taxon_id is returned, or None if no taxon is found. If *strict* is
            True and *rank* is not None, the set of taxon_ids found is
            resturned.

        Raises:
            KeyError: If multiple taxa are found, a *KeyError* exception is
            raised.
        """
        if not self._name_map:
            self.gen_name_map()

        if rank is None:
            return list(set(self._name_map[s_name.lower()]))
        else:
            taxon_ids = [
                taxon_id
                for taxon_id in self._name_map.get(s_name.lower(), [])
                if self[taxon_id].rank == rank
            ]
            if len(set(taxon_ids)) > 1:
                if strict:
                    raise KeyError(
                        "More than one taxon ({}) with rank {}".format(
                            s_name,
                            rank
                        )
                    )
                else:
                    return set(taxon_ids)
            elif len(taxon_ids) == 0:
                return None
            else:
                return taxon_ids[0]

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
        :return: instance of :class:`TaxonTuple` for the rank found.
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
            # kept only for backward compatibility.
            # needs a check to other scripts
            elif roots and (ranked.s_name in TAXON_ROOTS):
                break

            ranked = self[ranked.parent_id]

        return ranked

    def get_ranked_id(self, taxon_id, rank=None, it=False, include_higher=True):
        """
        .. versionadded:: 0.3.4

        Gets the ranked taxon of another one. Useful when it's better to get a
        taxon_id instead of an instance of :class:`TaxonTuple`. Internally, it
        relies on :meth:`Taxonomy.get_ranked_taxon`.

        Arguments:
            taxon_id (int): taxon_id
            rank (str or None): passed over
            it (bool): determines the return value. if True, a list is returned
            include_higher (bool): if True, any rank higher than the requested
                may be returned. If False and the rank cannot be returned, None
                is returned

        Returns:
            int or list: The type returned is based on the *it* paramenter. If
            *it* is True, the return value is a list with the *taxon_id* of the
            ranked taxon as the sole value. If False, the returned value is the
            *taxon_id*. *include_higher* determines if the return value should
            be *None* if the exact rank was not found and *include_higher* is
            False
        """
        taxon = self.get_ranked_taxon(taxon_id, rank)
        if (taxon.rank != rank) and (not include_higher):
            taxon = None
        else:
            taxon = taxon.taxon_id
        if it:
            taxon = [taxon]
        return taxon

    def get_name_map(self):
        "Returns a taxon_id->s_name dictionary"

        return dict(
            (taxon.taxon_id, taxon.s_name)
            for taxon in self
        )

    def get_lineage(self, taxon_id, names=False, only_ranked=True, with_last=True):
        """
        .. versionadded:: 0.3.1

        Proxy for :func:`get_lineage`, with changed defaults

        Arguments:
            taxon_id (int): taxon_id to return the lineage
            names (bool): if the elements of the list are converted into the
                scientific names
            only_ranked (bool): only return the ranked taxa
            with_last (bool): include the taxon_id passed to the list

        Returns:
            list: the lineage of the passed taxon_id as a list of IDs or names
        """

        return get_lineage(
            self,
            taxon_id,
            names=names,
            only_ranked=only_ranked,
            with_last=with_last
        )

    def get_lineage_string(self, taxon_id, only_ranked=True, with_last=True,
                           sep=';', rank=None):
        """
        .. versionadded:: 0.3.3

        Generates a lineage string, with the possibility of getting another
        ranked taxon (via :meth:`Taxonomy.get_ranked_taxon`) to another
        rank, such as *phylum*.

        Arguments:
            taxon_id (int): taxon_id to return the lineage
            only_ranked (bool): only return the ranked taxa
            with_last (bool): include the taxon_id passed to the list
            sep (str): separator used to join the lineage string
            rank (int or None): if None the full lineage is returned, otherwise
                the lineage will be cut to the specified rank

        Returns:
            str: lineage string
        """

        if rank is not None:
            taxon_id = self.get_ranked_taxon(taxon_id, rank=rank).taxon_id

        return sep.join(
            self.get_lineage(
                taxon_id,
                only_ranked=only_ranked,
                with_last=with_last,
                names=True
            )
        )

    def add_taxon(self, taxon_name, common_name='', rank='no rank', parent_id=None):
        """
        .. versionadded:: 0.3.1

        Adds a taxon to the taxonomy. If a taxon with the same name and rank is
        found, its taxon_id is returned, otherwise a new taxon_id is returned.

        Arguments:
            taxon_name (str): scientific name of the taxon
            common_name (str): common name
            rank (str): rank, defaults to 'no rank'
            parent_id (int): taxon_id of the parent, defaults to *None*, which
                is the taxonomy root

        Returns:
            int: the taxon_id of the added taxon (if new), or the taxon_id of
            the taxon with the same name and rank found in the taxonomy

        Raises:
            KeyError: if more than one taxon has already the passed name and
            rank and it can't be resolved by looking at the parent_id passed,
            the exception is raised.
        """
        # Tries to find it in the name map
        taxon_id = self.find_by_name(taxon_name, rank=rank, strict=False)

        # Multiple taxa with the same name, check if the parent_id of any of
        # them can be resolved
        if isinstance(taxon_id, set):
            # print taxon_id
            taxon_ids = taxon_id
            taxon_id = None
            for t_id in taxon_ids:
                if self[t_id].parent_id == parent_id:
                    taxon_id = t_id

            if taxon_id is None:
                raise KeyError(
                    "Multiple taxa ({}, {}) are in the taxonomy but cannot be associated with the parent_id ({}) provided".format(
                        taxon_name,
                        rank,
                        parent_id
                    )
                )

        if taxon_id is None:
            try:
                taxon_id = max(self._taxa) + 1
            except ValueError:
                taxon_id = 1

            self[taxon_id] = TaxonTuple(
                taxon_id,
                taxon_name,
                common_name,
                rank,
                (None,),
                parent_id
            )

            # adds the new name to the name map (which is initialised by
            # find_by_name)
            self._name_map[taxon_name.lower()] = [taxon_id]

        return taxon_id

    def add_lineage(self, **lineage):
        """
        .. versionadded:: 0.3.1

        Adds a lineage to the taxonomy. It's passed by keyword arguments, where
        each key is a value in the `TAXON_RANKS` rankes and the value is the
        scientific name. Appended underscores '_' will be stripped from the
        rank name. This is for cases such as *class* where the key is a reserved
        word in Python. Also one extra node can be added, such as
        strain/cultivar/subspecies and so on, but one only is expected to be
        passed.

        Arguments:
            lineage (dict): the lineage as a keyword arguments

        Returns:
            int: the taxon_id of the last element in the lineage

        Raises:
            ValueError: if more than a keyword argument is not contained in
            `TAXON_RANKS`
        """
        lineage = {
            key.rstrip('_'): value
            for key, value in viewitems(lineage)
        }
        extra_nodes = set(lineage) - set(TAXON_RANKS)
        if len(extra_nodes) > 1:
            raise ValueError(
                "Expected at most a one lower level below species, found: {}".format(
                    ', '.join(extra_nodes)
                )
            )
        else:
            extra_nodes = list(extra_nodes) if len(extra_nodes) == 1 else []

        ranks = list(TAXON_RANKS) + extra_nodes

        parent_ids = [None]
        taxon_id = None
        for rank in ranks:
            try:
                taxon_name = lineage[rank]
            except KeyError:
                continue

            if not isinstance(taxon_name, str):
                continue

            taxon_id = self.add_taxon(
                taxon_name,
                rank=rank,
                parent_id=parent_ids[-1]
            )
            # print taxon_id, rank, taxon_name, parent_ids[-1]
            parent_ids.append(taxon_id)

        return taxon_id

    def drop_taxon(self, taxon_id):
        """
        .. versionadded:: 0.3.1

        Drops a taxon and all taxa below it in the taxonomy. Also reset the
        name map for conistency.

        Arguments:
            taxon_id (int): taxon_id to drop from the taxonomy
        """
        drop_keys = []
        for leaf_id in self._taxa:
            if self.is_ancestor(leaf_id, taxon_id):
                drop_keys.append(leaf_id)

        for taxon_id in drop_keys:
            del self._taxa[taxon_id]

        self._name_map = {}

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
        for taxon in viewvalues(self._taxa):
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

    def __repr__(self):
        """
        .. versionadded:: 0.2.5
        """
        return "{} - {} taxa".format(self.__class__, len(self))


UniprotTaxonomy = Taxonomy


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

    :param :class:`Taxonomy` taxonomy: taxonomy used to test
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
        taxonomy: :class:`Taxonomy` instance used to test
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
        taxonomy: :class:`Taxonomy` instance
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
        taxonomy: :class:`Taxonomy` instance
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
        taxonomy: :class:`Taxonomy` instance
        taxon_ids (iterable): list taxonomic identifiers

    Returns:
        pandas.DataFrame: matrix with the pairwise distances of all *taxon_ids*
    """
    try:
        import pandas
    except ImportError:
        raise DependencyError('pandas')

    matrix = pandas.DataFrame(index=taxon_ids, columns=taxon_ids).fillna(0)

    for taxon_id1, taxon_id2 in itertools.combinations(taxon_ids, 2):
        distance = distance_two_taxa(taxonomy, taxon_id1, taxon_id2)
        matrix.loc[taxon_id1, taxon_id2] = distance
        matrix.loc[taxon_id2, taxon_id1] = distance

    return matrix


def get_lineage(taxonomy, taxon_id, names=False, only_ranked=False,
                with_last=False):
    """
    .. versionadded:: 0.2.1

    .. versionchanged:: 0.2.5
        added *only_ranked*

    .. versionchanged:: 0.3.0
        added *with_last*

    Returns the lineage of a taxon_id, as a list of taxon_id or taxa names

    Arguments:
        taxonomy: a :class:`Taxonomy` instance
        taxon_id (int): taxon_id whose lineage to return
        names (bool): if True, the returned list contains the names of the taxa
            instead of the taxon_id
        only_ranked (bool): if True, only taxonomic levels whose rank is in
            data:`TAXON_RANKS` will be returned
        with_last (bool): if True, the passed taxon_id is included in the
            lineage

    Returns:
        list: lineage of the taxon_id, the elements are `int` if names is False,
        and `str` when *names* is True. If a taxon has no scientific name, the
        common name is used. If *only_ranked* is True, the returned list only
        contains ranked taxa (according to :data:`TAXON_RANKS`).
    """
    lineage = []

    if with_last:
        lineage.append(taxon_id)

    while True:
        taxon_id = taxonomy[taxon_id].parent_id
        if taxon_id is None:
            break
        if only_ranked and (taxonomy[taxon_id].rank not in TAXON_RANKS):
            continue
        lineage.append(taxon_id)

    if names:
        lineage = [
            taxonomy[x].s_name if taxonomy[x].s_name else taxonomy[x].c_name
            for x in lineage
        ]

    return list(reversed(lineage))


def last_common_ancestor_multiple(taxonomy, taxon_ids):
    """
    .. versionadded:: 0.2.5

    Applies :func:`last_common_ancestor` to an iterable that yields *taxon_id*
    while removing any *None* values. If the list is of one element, that
    *taxon_id* is returned.

    Arguments:
        taxonomy: instance of :class:`Taxonomy`
        taxon_ids (iterable): an iterable that yields taxon_id

    Returns:
        int: the *taxon_id* that is the last common ancestor of all taxon_ids
        passed

    Raises:
        NoLcaFound: when no common ancestry is found or the number of
        *taxon_ids* is 0
    """

    taxon_ids = [
        taxon_id
        for taxon_id in taxon_ids
        if taxon_id is not None
    ]

    if not taxon_ids:
        raise NoLcaFound('Empty *taxon_ids*')

    return reduce(
        lambda x1, x2: last_common_ancestor(taxonomy, x1, x2),
        set(taxon_ids)
    )
