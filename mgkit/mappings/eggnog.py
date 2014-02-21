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

LOG = logging.getLogger(__name__)

EGGNOG_CAT_NAMES = (
    'INFORMATION STORAGE AND PROCESSING',
    'CELLULAR PROCESSES AND SIGNALING',
    'METABOLISM',
    'POORLY CHARACTERIZED'
)

EGGNOG_CAT_KEYS = (
    ('J', 'A', 'K', 'L', 'B'),
    ('D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'),
    ('C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q'),
    ('R', 'S')
)

EGGNOG_CAT_MAP = dict(
    (label, categories)
    for label, categories in zip(EGGNOG_CAT_NAMES, EGGNOG_CAT_KEYS)
)

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
        duplica funzionalita' in load_func_cat

        Le categorie sono una singola lettera e le desrizioni son nel dizionario EGGNOG_CAT,
        mentre le macro categorie si possono ricavare da EGGNOG_CAT_KEYS e EGGNOG_CAT_NAMES (con un ordine
        corrispondente)
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
