"""
Module containing classes and functions to deal with eggNOG data

.. todo::

    * unify download of data from web

"""
from __future__ import print_function

from .. import kegg
import urllib2
import pickle
import logging
import gzip
import cStringIO
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

    #INFORMATION STORAGE AND PROCESSING
    'J': "Translation, ribosomal structure and biogenesis",
    'A': "RNA processing and modification",
    'K': "Transcription",
    'L': "Replication, recombination and repair",
    'B': "Chromatin structure and dynamics",

    #CELLULAR PROCESSES AND SIGNALING
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

    #METABOLISM
    'C': "Energy production and conversion",
    'G': "Carbohydrate transport and metabolism",
    'E': "Amino acid transport and metabolism",
    'F': "Nucleotide transport and metabolism",
    'H': "Coenzyme transport and metabolism",
    'I': "Lipid transport and metabolism",
    'P': "Inorganic ion transport and metabolism",
    'Q': "Secondary metabolites biosynthesis, transport and catabolism",

    #POORLY CHARACTERIZED
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
        for gen_category, categories in EGGNOG_CAT_MAP.iteritems()
        if category in categories
    )


class Kegg2NogMapper(kegg.KeggMapperBase):
    """
    Usage

    to add a ko_id to the mapper use:
    * map_ko_to_eggnog(ko_id) which calls:

        * get_uniprot_from_ko
        * get_eggnog_from_uniprot

    return None if there was no problem or a dictionary with the exception thrown

    to make mappings ko_id -> eggnog-categories:
    * get_cat_mapping_from_file (file are in the current directory)
    * gen_ko_to_cat

    .. todo::

        * add method to load kos from file
        * optimise querying (like Kegg2CazyMapper)

    """
    query_string = 'database:(type:ko {0}) AND database:(type:eggnog)'
    columns_string = 'database(eggnog)'
    _egg_to_cat = None
    _ko_to_cat = None

    def __init__(self, fname=None, gen_ko_to_cat=True):
        super(Kegg2NogMapper, self).__init__()
        self._egg_to_cat = {}
        self._ko_to_cat = {}
        self._id_names = EGGNOG_CAT
        if fname:
            self.load_data(fname)
            if gen_ko_to_cat:
                self.gen_ko_to_cat()

    def map_ko_to_eggnog(self, ko_id, contact):
        if (ko_id in self._ko_map) or (ko_id in self._not_found):
            return
        try:
            mappings = self.ko_to_mapping(
                ko_id, self.query_string,
                self.columns_string, contact
            )
        except urllib2.HTTPError:
            return

        if mappings:
            self._ko_map[ko_id] = mappings
        else:
            self._not_found.append(ko_id)

    def map_kos_eggnog(self, kos, contact):
        errors = 0

        for idx, ko_id in enumerate(kos):
            LOG.info(
                "(%05d/%05d) Getting eggNOG mappings for KO: %s - %s", idx + 1,
                len(kos), ko_id, kos[ko_id]
            )
            try:
                self.map_ko_to_eggnog(ko_id, contact)
            except urllib2.HTTPError:
                LOG.warning("Problem getting mappings %s", ko_id)
                errors += 1

        return errors

    def get_cat_mapping_from_file(self, files=('COG.funccat.txt', 'NOG.funccat.txt', 'KOG.funccat.txt')):
        """
        Duplicate functionality of load_func_cat

        Categories are single letter and their descriptioons are in
        :data:`EGGNOG_CAT`, while the broader categories :data:`EGGNOG_CAT_MAP`
        """
        for fname in files:
            f = open(fname, 'r')
            for line in f:
                egg_id, cats = line.rstrip().split()
                self._egg_to_cat[egg_id] = list(cats)

    def get_cat_mapping(self, f_handle):
        for line in f_handle:
            line = line.rstrip()
            if line == '':
                continue
            egg_id, cats = line.rstrip().split('\t', 1)
            self._egg_to_cat[egg_id] = list(cats)

    def gen_ko_to_cat(self):
        for ko_id, egg_ids in self._ko_map.iteritems():
            cats = set()
            for egg_id in egg_ids:
                try:
                    cats.update(self._egg_to_cat[egg_id])
                except KeyError:
                    #some egg_ids have no assigned category
                    pass
            if cats:
                self._ko_to_cat[ko_id] = cats

    def get_ko_cat(self, ko_id):
        return sorted(self._ko_to_cat[ko_id])

    def get_ko_cat_dict(self):
        return self._ko_to_cat.copy()

    def get_ko_map(self):
        return self.get_ko_cat_dict()

    def get_cat_ko_dict(self):
        cat_map = {}

        for ko_id in self._ko_to_cat:
            for mapped_id in self._ko_to_cat[ko_id]:
                try:
                    cat_map[mapped_id].add(ko_id)
                except KeyError:
                    cat_map[mapped_id] = set([ko_id])
        return cat_map

    def save_data(self, fname):
        LOG.info("Saving data to %s", fname)
        pickle.dump(
            (self._ko_map, self._egg_to_cat, self._not_found), open(fname, 'w')
        )

    def load_data(self, fname):
        LOG.info("Loading data from %s", fname)
        self._ko_map, self._egg_to_cat, self._not_found = pickle.load(
            open(fname, 'r')
        )


def download_data(contact, kegg_data='kegg.pickle', eggnog_data='eggnog.pickle',
                  base_url='http://eggnog.embl.de/version_3.0/data/downloads/'):

    egg = Kegg2NogMapper()
    try:
        egg.load_data(eggnog_data)
    except:
        LOG.warning("Eggnog data not found, download restarting")

    kdata = kegg.KeggData(kegg_data)
    ko_names = kdata.get_ko_names()

    urls = ('COG.funccat.txt.gz', 'KOG.funccat.txt.gz', 'NOG.funccat.txt.gz')

    try:
        errors = egg.map_kos_eggnog(ko_names, contact)

        for url in urls:
            LOG.info("Loading categories mapping from %s%s", base_url, url)
            f_handle = gzip.GzipFile(
                fileobj=cStringIO.StringIO(urllib2.urlopen(base_url + url).read()),
                mode='r'
            )
            egg.get_cat_mapping(f_handle)
    finally:
        egg.save_data(eggnog_data)

    LOG.info("Number of errors: %d", errors)


def print_ko_to_cat(data_file='eggnog.pickle', descr=False):
    e = Kegg2NogMapper()
    e.load_data(data_file)
    e.gen_ko_to_cat()

    for ko_id in sorted(e):
        try:
            cats = e.get_ko_cat(ko_id)
        except KeyError:
            cats = []

        print("{0}\t{1}".format(ko_id, '\t'.join(
            EGGNOG_CAT[x] if descr else x
            for x in cats)
        ))


class NOGInfo(object):
    """
    .. versionadded:: 0.1.14

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

    def __init__(self):
        self._map_nog_gene = {}
        self._map_gene_nog = {}
        self._map_nog_desc = {}
        self._map_nog_func = {}

    def load_members(self, file_handle):
        """
        Loads data from *NOG.members.txt.gz*

        *file_handle* can either an open file or a path
        """
        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'r')

        LOG.info("Reading Members from %s", file_handle.name)

        map_nog_gene = {}

        for line in file_handle:

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
        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'r')

        LOG.info("Reading Descriptions from %s", file_handle.name)

        map_nog_desc = {}

        for line in file_handle:

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
        if isinstance(file_handle, str):
            file_handle = open_file(file_handle, 'r')

        LOG.info("Reading Functional Categories from %s", file_handle.name)

        map_nog_func = {}

        for line in file_handle:
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
        for the requested eggNOG COG/NOG ID
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
