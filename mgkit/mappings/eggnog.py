"""
Module containing classes and functions to deal with eggNOG data

.. todo::

    * unify download of data from web

"""
from __future__ import print_function
from future.utils import viewitems
import logging
import itertools
from ..io import open_file
from ..utils import dictionary

LOG = logging.getLogger(__name__)

EGGNOG_CAT_NAMES = (
    'INFORMATION STORAGE AND PROCESSING',
    'CELLULAR PROCESSES AND SIGNALING',
    'METABOLISM',
    'POORLY CHARACTERIZED'
)
"""
Functional categories (broader)
"""

EGGNOG_CAT_KEYS = (
    ('J', 'A', 'K', 'L', 'B'),
    ('D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'),
    ('C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'),
    ('R', 'S')
)
"""
Used to build map of broader categories (:data:`EGGNOG_CAT_NAMES`) to more
specific ones
"""

EGGNOG_CAT_MAP = dict(
    (label, categories)
    for label, categories in zip(EGGNOG_CAT_NAMES, EGGNOG_CAT_KEYS)
)
"""
Functional categories (broader, :data:`EGGNOG_CAT_NAMES`) mappings to more
specific one (:data:`EGGNOG_CAT`).
"""

EGGNOG_CAT = {

    # INFORMATION STORAGE AND PROCESSING
    'J': "Translation, ribosomal structure and biogenesis",
    'A': "RNA processing and modification",
    'K': "Transcription",
    'L': "Replication, recombination and repair",
    'B': "Chromatin structure and dynamics",

    # CELLULAR PROCESSES AND SIGNALING
    'D': "Cell cycle control, cell division, chromosome partitioning",
    'Y': "Nuclear structure",
    'V': "Defense mechanisms",
    'T': "Signal transduction mechanisms",
    'M': "Cell wall/membrane/envelope biogenesis",
    'N': "Cell motility",
    'Z': "Cytoskeleton",
    'W': "Extracellular structures",
    'U': "Intracellular trafficking, secretion, and vesicular transport",
    'O': "Posttranslational modification, protein turnover, chaperones",

    # METABOLISM
    'C': "Energy production and conversion",
    'G': "Carbohydrate transport and metabolism",
    'E': "Amino acid transport and metabolism",
    'F': "Nucleotide transport and metabolism",
    'H': "Coenzyme transport and metabolism",
    'I': "Lipid transport and metabolism",
    'P': "Inorganic ion transport and metabolism",
    'Q': "Secondary metabolites biosynthesis, transport and catabolism",

    # POORLY CHARACTERIZED
    'R': "General function prediction only",
    'S': "Function unknown"
}
"""
Single letter functional categories
"""


def get_general_eggnog_cat(category):
    """
    .. versionadded:: 0.1.14

    Returns the functional category (:data:`EGGNOG_CAT_NAMES` keys)
    for the requested single letter functional category (:data:`EGGNOG_CAT`
    keys)
    """
    return set(
        gen_category
        for gen_category, categories in viewitems(EGGNOG_CAT_MAP)
        if category in categories
    )


class NOGInfo(object):
    """
    .. versionadded:: 0.1.14

    .. versionchanged:: 0.4.0
        made file reading compatible with Python 3

    Mappings from Uniprot to eggNOG

    ..note::

        load_description is optional
    """
    _map_nog_func = None
    "eggNOG COG/NOG to functional category dictionary"
    _map_nog_gene = None
    "eggNOG COG/NOG to gene id dictionary"
    _map_nog_desc = None
    "eggNOG id to description dictionary"
    _map_gene_nog = None
    "eggNOG gene id to COG/NOG dictionary"

    def __init__(self, members=None, funccat=None, description=None):
        """
        .. versionchanged:: 0.4.0
            added parameters to load data at __init__
        """
        self._map_nog_gene = {}
        self._map_gene_nog = {}
        self._map_nog_func = {}
        self._map_nog_desc = {}
        if members is not None:
            self.load_members(members)
        if funccat is not None:
            self.load_funccat(funccat)
        if description is not None:
            self.load_description(description)

    def load_members(self, file_handle):
        """
        Loads data from *NOG.members.txt.gz*

        *file_handle* can either an open file or a path
        """
        file_handle = open_file(file_handle, 'rb')

        LOG.info("Reading Members from %s", file_handle.name)

        map_nog_gene = {}

        for line in file_handle:
            line = line.decode('ascii')
            if line.startswith('#'):
                continue

            nog_id, gene_id, start, end = line.strip().split()

            try:
                map_nog_gene[nog_id].append(gene_id)
            except KeyError:
                map_nog_gene[nog_id] = [gene_id]

        self._map_nog_gene.update(map_nog_gene)
        self._map_gene_nog.update(dictionary.reverse_mapping(map_nog_gene))

    def load_description(self, file_handle):
        """
        Loads data from *NOG.description.txt.gz*

        *file_handle* can either an open file or a path
        """
        file_handle = open_file(file_handle, 'rb')

        LOG.info("Reading Descriptions from %s", file_handle.name)

        map_nog_desc = {}

        for line in file_handle:
            line = line.decode('ascii')
            try:
                nog_id, nog_desc = line.strip().split('\t', 1)
            except ValueError:
                nog_id = line.strip()
                nog_desc = ''

            map_nog_desc[nog_id] = nog_desc

        self._map_nog_desc.update(map_nog_desc)

    def load_funccat(self, file_handle):
        """
        Loads data from *NOG.funccat.txt.gz*

        *file_handle* can either an open file or a path
        """
        file_handle = open_file(file_handle, 'rb')

        LOG.info("Reading Functional Categories from %s", file_handle.name)

        map_nog_func = {}

        for line in file_handle:
            line = line.decode('ascii')
            nog_id, nog_func = line.strip().split()
            nog_func = set(nog_func)

            map_nog_func[nog_id] = nog_func

        self._map_nog_func.update(map_nog_func)

    def get_nog_funccat(self, nog_id):
        """
        Returns the functional category (one letter, :data:`EGGNOG_CAT` keys)
        for the requested eggNOG COG/NOG ID
        """
        try:
            return self._map_nog_func[nog_id].copy()
        except KeyError:
            return set()

    def get_nogs_funccat(self, nog_ids):
        """
        Returns the functional categories for a list of COG/NOG IDs. Uses
        :meth:`NOGInfo.get_nog_funccat`
        """

        iterator = (self.get_nog_funccat(nog_id) for nog_id in nog_ids)

        return set(itertools.chain(*iterator))

    def get_nog_gencat(self, nog_id):
        """
        Returns the functional category (:data:`EGGNOG_CAT_NAMES` keys)
        for the requested eggNOG COG/NOG IDs
        """
        return set(
            itertools.chain(
                *(
                    get_general_eggnog_cat(funccat)
                    for funccat in self.get_nog_funccat(nog_id)
                )
            )
        )

    def get_gene_nog(self, gene_id):
        """
        Returns the COG/NOG ID of the requested eggNOG gene ID
        """
        try:
            nog_ids = self._map_gene_nog[gene_id]
        except KeyError:
            return None

        return nog_ids

    def get_gene_funccat(self, gene_id):
        """
        Returns the functional category (one letter, :data:`EGGNOG_CAT` keys)
        for the requested eggNOG gene ID
        """
        try:
            nog_ids = self._map_gene_nog[gene_id]
        except KeyError:
            return None

        return self.get_nogs_funccat(nog_ids)
