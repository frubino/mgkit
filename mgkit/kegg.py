"""
Module containing classes and functions to access Kegg data
"""

import logging
import pickle
import re
import itertools
from .utils import dictionary as dict_utils
from .utils.common import deprecated
from .net import uniprot, url_read


LOG = logging.getLogger(__name__)
KEGG_REST_URL = 'http://rest.kegg.jp/'


class KeggCompound(object):
    "Kegg compound"
    # __slots__ = ('cp_id', 'description')

    def __init__(self, cp_id=None, description=''):
        self.cp_id = cp_id
        self.description = description

    def __eq__(self, other):
        """
        >>> KeggCompound('test') == KeggCompound('test')
        True
        >>> KeggCompound('test') == 1
        False
        """
        if not isinstance(other, KeggCompound):
            return False
        return self.cp_id == other.cp_id

    def __ne__(self, other):
        """
        >>> KeggCompound('test') != KeggCompound('test1')
        True
        >>> KeggCompound('test') != 1
        True
        """
        return not self == other

    def __hash__(self):
        return hash(self.cp_id)

    def __str__(self):
        return "{0}: {1}".format(self.cp_id, self.description)

    def __repr__(self):
        return str(self)


class KeggReaction(object):
    "Kegg Reaction"
    # __slots__ = ('rn_id', 'description', 'cp_in', 'cp_out')

    def __init__(self, rn_id=None, description='', cp_in=None, cp_out=None):
        self.rn_id = rn_id
        self.description = description
        self.cp_in = cp_in if not cp_in is None else {}
        self.cp_out = cp_out if not cp_out is None else {}

    def __eq__(self, other):
        """
        >>> KeggReaction('test') == KeggReaction('test')
        True
        >>> KeggReaction('test') == 1
        False
        """
        if not isinstance(other, KeggReaction):
            return False
        return self.rn_id == other.rn_id

    def __ne__(self, other):
        """
        >>> KeggReaction('test') != KeggReaction('test1')
        True
        >>> KeggReaction('test') != 1
        True
        """
        return not self == other

    def __hash__(self):
        return hash(self.rn_id)

    def __str__(self):
        return "{0}: {1} - CPIN ({2}) CPOUT ({3})".format(
            self.rn_id, self.description, len(self.cp_in), len(self.cp_out)
        )

    def __repr__(self):
        return str(self)


class KeggOrtholog(object):
    "Kegg Ortholog gene"
    #__slots__ = ('ko_id', 'description', 'reactions')

    def __init__(self, ko_id=None, description='', reactions=None):
        self.ko_id = ko_id
        self.description = description.replace(';', '')
        self.reactions = reactions if not reactions is None else {}

    def __eq__(self, other):
        """
        >>> KeggOrtholog('test') == KeggOrtholog('test')
        True
        >>> KeggOrtholog('test') == 1
        False
        """
        if not isinstance(other, KeggOrtholog):
            return False
        return self.ko_id == other.ko_id

    def __ne__(self, other):
        """
        >>> KeggOrtholog('test') != KeggOrtholog('test1')
        True
        >>> KeggOrtholog('test') != 1
        True
        """
        return not self == other

    def __getitem__(self, key):
        return self.reactions[key]

    def __setitem__(self, key, value):
        self.reactions[key] = value

    def __hash__(self):
        return hash(self.ko_id)

    def __len__(self):
        return len(self.reactions)

    def __iter__(self):
        for x in self.reactions:
            yield x

    def __str__(self):
        return "{0}: ({1}) {2}".format(self.ko_id, len(self), self.description)

    def __repr__(self):
        return str(self)


class KeggPathway(object):
    "Kegg Pathway"
    # __slots__ = ('path_id', 'description', 'genes')

    def __init__(self, path_id=None, description=None, genes=None):
        self.path_id = path_id
        self.description = description
        self.genes = genes if not genes is None else {}

    def __eq__(self, other):
        """
        >>> KeggPathway('test') == KeggPathway('test')
        True
        >>> KeggPathway('test') == 1
        False
        """
        if not isinstance(other, KeggPathway):
            return False
        return self.path_id == other.path_id

    def __ne__(self, other):
        """
        >>> KeggPathway('test') != KeggPathway('test1')
        True
        >>> KeggPathway('test') != 1
        True
        """
        return not self == other

    def __hash__(self):
        return hash(self.path_id)

    def __contains__(self, item):
        return item in self.genes

    def __getitem__(self, key):
        return self.genes[key]

    def __setitem__(self, key, value):
        self.genes[key] = value

    def __str__(self):
        return "{0} ({2}): {1}".format(
            self.path_id, self.description, len(self)
        )

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.genes)

    def __iter__(self):
        for x in self.genes:
            yield x


class KeggClientRest(object):
    """
    Kegg REST client

    The class includes methods and data to use the REST API provided by Kegg. At
    the moment it provides methods to for 'link', 'list' and 'get' operations,

    `Kegg REST API <http://www.kegg.jp/kegg/rest/keggapi.html>`_

    """
    contact = None
    api_url = KEGG_REST_URL
    cpd_re = re.compile(
        r"ENTRY\s+(C\d{5})\s+Compound\nNAME\s+([,.\w+ ()-]+);?"
    )
    rn_name_re = re.compile(r"R\d{5}")
    rn_eq_re = re.compile(r'C\d{5}')
    ko_desc_re = re.compile(
        r"ko:(K\d{5})\t.+?;\s+([\w+, ()/:'\[\]-]+)( \[EC:)?\n?"
    )
    cpd_desc_re = re.compile(
        r"cpd:(C\d{5})\t([\w+, ()\[\]'*.-]+);?\n?"
    )
    id_prefix = {'C': 'cpd', 'k': 'map', 'K': 'ko', 'R': 'rn', 'm': 'path'}

    ####### Kegg primitives #######

    def link_ids(self, target, ids, strip=True, max_len=50):
        """
        The method abstract the use of the 'link' operation in the Kegg API

        The target parameter can be one of the following::

            pathway | brite | module | disease | drug | environ | ko | genome |
            <org> | compound | glycan | reaction | rpair | rclass | enzyme

            <org> = KEGG organism code or T number

        :param str target: the target db
        :param ids: can be either a single id as a string or a list of ids
        :param bool strip: if the prefix (e.g. ko:K00601) should be stripped
        :param int max_len: the maximum number of ids to retrieve with each
            request, should not exceed 50
        :return dict: dictionary mapping requested id to target id(s)
        """
        if isinstance(ids, str):
            ids = [ids]
        if isinstance(ids, (set, frozenset, dict)):
            ids = list(ids)

        data = []
        for idx in range(0, len(ids), max_len):
            if len(ids) > max_len:
                LOG.info(
                    "Downloading links - range %d-%d",
                    idx + 1, idx + max_len
                )
            url = "{0}link/{1}/{2}".format(
                self.api_url, target, '+'.join(ids[idx:idx+max_len])
            )
            t_data = url_read(url, agent=self.contact)
            t_data = t_data.split('\n')
            data.extend(t_data)

        mapping = {}
        for line in data:
            if not line:
                continue
            source, target = line.split()
            if strip:
                source = source.split(':')[1]
                target = target.split(':')[1]
            try:
                mapping[source].append(target)
            except KeyError:
                mapping[source] = [target]

        return mapping

    def list_ids(self, k_id):
        """
        The method abstract the use of the 'list' operation in the Kegg API

        The k_id parameter can be one of the following::

            pathway | brite | module | disease | drug | environ | ko | genome |
            <org> | compound | glycan | reaction | rpair | rclass | enzyme

            <org> = KEGG organism code or T number

        :param str k_id: kegg database to get list of ids
        :return list: list of ids in the specified database
        """
        url = "{0}list/{1}".format(self.api_url, k_id)
        data = url_read(url, agent=self.contact)
        #leave out the last \n
        return data[:-1]

    def get_entry(self, k_id, option=None):
        """
        The method abstract the use of the 'get' operation in the Kegg API

        :param str k_id: kegg id of the resource to get
        :param str option: optional, to specify a format
        """
        url = "{0}get/{1}/{2}".format(
            self.api_url,
            k_id,
            '' if option is None else option
        )
        data = url_read(url, agent=self.contact)
        return data

    def link(self, target, source, options=None):
        """
        .. versionadded:: 0.2.0

        Implements "link" operation in Kegg REST

        http://www.genome.jp/linkdb/
        """
        url = "http://rest.genome.jp/link/{0}/{1}/{2}".format(
            target,
            source,
            '' if options is None else options
        )

        LOG.debug(url)

        data = url_read(url, agent=self.contact)

        mappings = {}
        for line in data.splitlines():
            source_id, target_id, _ = line.rstrip().split('\t')
            source_id = source_id.split(':')[1]
            target_id = target_id.split(':')[1]
            try:
                mappings[source_id].add(target_id)
            except KeyError:
                mappings[source_id] = set([target_id])

        return mappings

    ####### names #######

    def get_ids_names(self, target='ko', strip=True):
        """
        .. versionadded:: 0.1.13

        Returns a dictionary with the names/description of all the id of a
        specific target, (ko, path, cpd, etc.)

        If strip=True the id will stripped of the module abbreviation (e.g.
        md:M00002->M00002)
        """
        id_names = {}

        for line in self.list_ids(target).splitlines():
            kegg_id, name = line.strip().split('\t')
            if strip:
                kegg_id = kegg_id.split(':')[1]
            id_names[kegg_id] = name
        return id_names

    @deprecated
    def get_names(self, k_type, strip=True):
        """
        .. deprecated:: 0.1.13
            use :meth:`KeggClientRest.get_ids_names`

        Returns a dictionary kegg_id->description for the specified kegg module.

        If strip=True the id will stripped of the module abbreviation (e.g.
        md:M00002->M00002)
        """
        return self.get_ids_names(target=k_type, strip=strip)

    @deprecated
    def get_compound_names(self, ids=None):
        """
        .. deprecated:: 0.1.13
            use :meth:`KeggClientRest.get_ids_names`
        """
        if ids is not None:
            return self.get_compound_names(ids)

        return self.get_ids_names('cpd')

    @deprecated
    def get_kos_descriptions(self, rex=False):
        """
        .. deprecated:: 0.1.13
            use :meth:`KeggClientRest.get_ids_names`

        Get all KOs descriptions

        """

        return self.get_ids_names(target='ko')

    def get_ortholog_pathways(self):
        """
        Gets ortholog pathways, replace 'map' with 'ko' in the id
        """
        data = self.get_ids_names('pathway')
        pathways = {}
        for kegg_id, name in data.iteritems():

            kegg_id = kegg_id.replace('map', 'ko')
            pathways[kegg_id] = name

        return pathways

    @deprecated
    def get_compounds_descriptions(self):
        """
        .. deprecated:: 0.1.13
            use *get_compound_names* instead

        Get compound descriptions
        """
        return self.get_ids_names('cpd')

    @deprecated
    def get_reactions_descriptions(self):
        """
        .. deprecated:: 0.1.13
            use *get_compound_names* instead

        Get reaction descriptions
        """
        return self.get_ids_names('reaction')

    ####### end names #######

    def get_reaction_equations(self, ids, max_len=10):
        "Get the equation for the reactions"
        if isinstance(ids, str):
            ids = [ids]
        data = {}
        for idx in range(0, len(ids), max_len):
            if len(ids) > max_len:
                LOG.info(
                    "Downloading reactions equations - range %d-%d",
                    idx + 1, idx + max_len)

            url = "{0}get/{1}".format(
                self.api_url, '+'.join(ids[idx:idx+max_len])
            )

            t_data = url_read(url)
            for entry in t_data.split('///'):
                for line in entry.split('\n'):
                    if line.startswith('EQUATION'):
                        cp_in, cp_out = line.split('<=>')
                        # print cp_in, cp_out
                        cp_in = self.rn_eq_re.findall(cp_in)
                        cp_out = self.rn_eq_re.findall(cp_out)
                        # print cp_in, cp_out
                    elif line.startswith('ENTRY'):
                        # try:
                        name = self.rn_name_re.search(line).group(0)
                        # except:
                        #     print line
                data[name] = {'in': cp_in, 'out': cp_out}

        return data

    def get_pathway_links(self, pathway):
        """
        Returns a dictionary with the mappings KO->compounds for a specific
        Pathway or module
        """
        compounds = self.link_ids('cpd', pathway)
        compounds = set(compounds[pathway])
        kos = self.link_ids('ko', pathway)
        kos = set(kos[pathway])
        cp_rn = self.link_ids('ec', list(compounds))
        ko_rn = self.link_ids('ec', list(kos))
        cp_rn_rev = dict_utils.reverse_mapping(cp_rn)
        ko_cp = dict_utils.combine_dict(ko_rn, cp_rn_rev)

        return ko_cp


class KeggData(object):
    pathways = None
    _ko_map = None
    maps = None

    def __init__(self, fname=None, gen_maps=True):
        self.pathways = {}
        self._ko_map = {}
        self.maps = None
        if fname:
            self.load_data(fname)
            if gen_maps:
                self.gen_ko_map()
                self.gen_maps()

    def get_ko_names(self):
        ko_names = {}
        for path_id in self:
            for ko_id in self[path_id]:
                if ko_id in ko_names:
                    continue
                ko_names[ko_id] = self[path_id][ko_id].description
        return ko_names

    def get_rn_names(self):
        rn_names = {}
        for path_id in self:
            for ko_id in self[path_id]:
                for rn_id in self[path_id][ko_id].reactions:
                    if rn_id in rn_names:
                        continue
                    rn_names[rn_id] = self[path_id][ko_id].reactions[rn_id].description
        return rn_names

    def get_cp_names(self):
        if self.maps is None:
            self.gen_maps()
        return dict(
            (cp_id, cp.description)
            for cp_id, cp in self.maps['cp'].iteritems()
        )

    def get_ko_rn_links(self, path_filter=None, description=False):
        rn_links = {}
        for path_id in self:
            if not path_filter is None:
                if path_id != path_filter:
                    continue
            for ko_id in self[path_id]:
                if description:
                    rn_links[ko_id] = [
                        "{0}: {1}".format(rn.rn_id, rn.description)
                        for rn in self[path_id][ko_id].reactions.itervalues()
                    ]
                else:
                    rn_links[ko_id] = self[path_id][ko_id].reactions.keys()
        return rn_links

    def get_rn_cp_links(self, path_filter=None, description=False):
        rn_links = {}
        for path_id in self:
            if not path_filter is None:
                if path_id != path_filter:
                    continue
            for ko_id in self[path_id]:
                for rn in self[path_id][ko_id].reactions.itervalues():
                    if rn.rn_id in rn_links:
                        continue
                    cp_set = set()
                    if description:
                        cp_set.update(
                            "{0}: {1}".format(cp.cp_id, cp.description)
                            for cp in itertools.chain(
                                rn.cp_in.itervalues(), rn.cp_out.itervalues()
                            )
                        )
                    else:
                        cp_set.update(rn.cp_in, rn.cp_out)
                    rn_links[rn.rn_id] = cp_set
        return rn_links

    def get_ko_cp_links(self, path_filter=None, description=False):
        ko_links = {}
        for path_id in self:
            if not path_filter is None:
                if path_id != path_filter:
                    continue
            for ko_id in self[path_id]:
                if ko_id in ko_links:
                    continue
                cp_set = set()
                for rn in self[path_id][ko_id].reactions.itervalues():
                    if description:
                        cp_set.update(
                            "{0}: {1}".format(
                            cp.cp_id, cp.description)
                            for cp in itertools.chain(
                                rn.cp_in.itervalues(), rn.cp_out.itervalues()
                            )
                        )
                    else:
                        cp_set.update(rn.cp_in, rn.cp_out)
                ko_links[ko_id] = cp_set
        return ko_links

    def get_ko_cp_links_alt(self, direction='out', description=False):
        cp_links = {}
        for path_id in self:
            for ko_id in self[path_id]:
                if ko_id in cp_links:
                    continue
                cp = set()
                for rn in self[path_id][ko_id].reactions.itervalues():
                    cp_dict = getattr(rn, 'cp_' + direction)
                    if description:
                        for cp_obj in cp_dict.itervalues():
                            cp.add("{0}: {1}".format(
                                   cp_obj.cp_id, cp_obj.description))
                    else:
                        cp.update(cp_dict.keys())
                cp_links[ko_id] = cp
        return cp_links

    def gen_ko_map(self):
        for path_id in self:
            for ko_id in self[path_id]:
                try:
                    self._ko_map[ko_id].add(path_id)
                except KeyError:
                    self._ko_map[ko_id] = set()
                    self._ko_map[ko_id].add(path_id)

    def gen_maps(self):
        ko_maps = {}
        rn_maps = {}
        cp_maps = {}
        for path_id in self:
            for ko in self[path_id].genes.itervalues():
                if ko.ko_id in ko_maps:
                    continue
                ko_maps[ko.ko_id] = ko
                for rn in ko.reactions.itervalues():
                    if rn.rn_id in rn_maps:
                        continue
                    rn_maps[rn.rn_id] = rn
                    for cp in itertools.chain(rn.cp_in.values(),
                                              rn.cp_out.values()):
                        if cp.cp_id in cp_maps:
                            continue
                        cp_maps[cp.cp_id] = cp
        self.maps = {'cp': cp_maps, 'ko': ko_maps, 'rn': rn_maps}

    def get_pathway_ko_map(self, black_list=None):
        if black_list is None:
            black_list = self.pathways.keys()
        else:
            #keeps only path_ids that are not in the black_list
            black_list = set(self.pathways.keys()) - set(black_list)

        return dict(
            (path_id, self[path_id].genes.keys())
            for path_id in self
            if path_id in black_list
        )

    def get_ko_pathway_map(self, black_list=None):
        return dict_utils.reverse_mapping(
            self.get_pathway_ko_map(black_list=black_list)
        )

    def get_ko_pathways(self, ko_id):
        return sorted(self._ko_map[ko_id])

    def save_data(self, fname):
        LOG.info("Saving data to file %s", fname)
        pickle.dump(self.pathways, open(fname, 'w'))

    def load_data(self, fname):
        LOG.info("Loading data from file %s", fname)
        self.pathways = pickle.load(open(fname, 'r'))

    def __getitem__(self, key):
        return self.pathways[key]

    def __len__(self):
        return len(self.pathways)

    def __iter__(self):
        for x in self.pathways:
            yield x


class KeggMapperBase(object):
    """
    Base object for Kegg mapping classes
    """
    _ko_map = None
    _not_found = None
    _id_names = None

    def __init__(self, fname=None):
        self._ko_map = {}
        self._not_found = []
        self._id_names = {}
        if fname:
            self.load_data(fname)

    ko_to_mapping = staticmethod(uniprot.ko_to_mapping)

    def save_data(self, fname):
        """
        Saves mapping data to disk
        """
        LOG.info("Saving data to %s", fname)
        pickle.dump((self._ko_map, self._not_found, self._id_names),
                    open(fname, 'w'))

    def load_data(self, fname):
        """
        Loads mapping data to disk
        """
        LOG.info("Loading data from %s", fname)
        data = pickle.load(open(fname, 'r'))
        try:
            self._ko_map, self._not_found, self._id_names = data
        except ValueError:
            LOG.warning("Pickled data in old format")
            self._ko_map, self._not_found = data

    def get_ko_map(self):
        """
        Returns a copy of the KO->mapping dictionary
        """
        return self._ko_map.copy()

    def get_id_map(self):
        """
        Returns a mapping->KOs dictionary (a reverse mapping to get_ko_map)
        """
        return dict_utils.reverse_mapping(self.get_ko_map())

    def get_id_names(self):
        """
        Returns a copy of the mapping names
        """
        return self._id_names.copy()

    def __len__(self):
        return len(self._ko_map)

    def __getitem__(self, key):
        return self._ko_map[key]

    def __setitem__(self, key, value):
        self._ko_map[key] = value

    def __delitem__(self, key):
        del self._ko_map[key]

    def __contains__(self, key):
        return key in self._ko_map

    def __iter__(self):
        return iter(self._ko_map)

BLACK_LIST = [
    'ko05164',  # Influenza A
    'ko05166',  # HTLV-I infection
    'ko05161',  # Hepatitis B
    'ko05160',  # Hepatitis C
    'ko05162',  # Measles
    'ko05169',  # Epstein-Barr virus infection
    'ko05168',  # Herpes simplex infection
    'ko04520',  # Adherens junction
    'ko05014',  # Amyotrophic lateral sclerosis (ALS)
    'ko05016',  # Huntington's disease
    'ko05310',  # Asthma
    'ko00908',  # Zeatin biosynthesis
    'ko04930',  # Type II diabetes mellitus
    'ko04622',  # RIG-I-like receptor signaling pathway
    'ko04626',  # Plant-pathogen interaction
    'ko05222',  # Small cell lung cancer
    'ko05223',  # Non-small cell lung cancer
    'ko05220',  # Chronic myeloid leukemia
    'ko05221',  # Acute myeloid leukemia
    'ko04612',  # Antigen processing and presentation
    'ko04621',  # NOD-like receptor signaling pathway
    'ko04620',  # Toll-like receptor signaling pathway
    'ko04623',  # Cytosolic DNA-sensing pathway
    'ko04622',  # RIG-I-like receptor signaling pathway
    'ko00351',  # DDT degradation
    'ko04270',  # Vascular smooth muscle contraction
    'ko04115',  # p53 signaling pathway
    'ko04114',  # Oocyte meiosis
    'ko05120',  # Epithelial cell signaling in Helicobacter pylori infection
    'ko00590',  # Arachidonic acid metabolism
    'ko04370',  # VEGF signaling pathway
    'ko04976',  # Bile secretion
    'ko04971',  # Gastric acid secretion
    'ko04970',  # Salivary secretion
    'ko04973',  # Carbohydrate digestion and absorption
    'ko04972',  # Pancreatic secretion
    'ko04975',  # Fat digestion and absorption
    'ko04974',  # Protein digestion and absorption
    'ko04977',  # Vitamin digestion and absorption
    'ko04711',  # Circadian rhythm - fly
    'ko04710',  # Circadian rhythm
    'ko04713',  # Circadian entrainment
    'ko04712',  # Circadian rhythm - plant
    'ko04650',  # Natural killer cell mediated cytotoxicity
    'ko04391',  # Hippo signaling pathway - fly
    'ko04151',  # PI3K-Akt signaling pathway
    'ko04150',  # mTOR signaling pathway
    'ko04510',  # Focal adhesion
    'ko04012',  # ErbB signaling pathway
    'ko04013',  # MAPK signaling pathway - fly
    'ko04010',  # MAPK signaling pathway
    'ko03320',  # PPAR signaling pathway
    'ko04530',  # Tight junction
    'ko00981',  # Insect hormone biosynthesis
    'ko00982',  # Drug metabolism - cytochrome P450
    'ko00983',  # Drug metabolism - other enzymes
    'ko05150',  # Staphylococcus aureus infection
    'ko05152',  # Tuberculosis
    'ko04614',  # Renin-angiotensin system
    'ko04610',  # Complement and coagulation cascades
    'ko04340',  # Hedgehog signaling pathway
    'ko00984',  # Steroid degradation
    'ko05010',  # Alzheimer's disease
    'ko05012',  # Parkinson's disease
    'ko04260',  # Cardiac muscle contraction
    'ko04540',  # Gap junction
    'ko05110',  # Vibrio cholerae infection
    'ko05111',  # Vibrio cholerae pathogenic cycle
    'ko04640',  # Hematopoietic cell lineage
    'ko04310',  # Wnt signaling pathway
    'ko04728',  # Dopaminergic synapse
    'ko04724',  # Glutamatergic synapse
    'ko04725',  # Cholinergic synapse
    'ko04726',  # Serotonergic synapse
    'ko04727',  # GABAergic synapse
    'ko04720',  # Long-term potentiation
    'ko04721',  # Synaptic vesicle cycle
    'ko04722',  # Neurotrophin signaling pathway
    'ko04723',  # Retrograde endocannabinoid signaling
    'ko05340',  # Primary immunodeficiency
    'ko05414',  # Dilated cardiomyopathy
    'ko05416',  # Viral myocarditis
    'ko05410',  # Hypertrophic cardiomyopathy (HCM)
    'ko05412',  # Arrhythmogenic right ventricular cardiomyopathy (ARVC)
    'ko04130',  # SNARE interactions in vesicular transport
    'ko05145',  # Toxoplasmosis
    'ko05144',  # Malaria
    'ko04350',  # TGF-beta signaling pathway
    'ko04210',  # Apoptosis
    'ko05200',  # Pathways in cancer
    'ko05202',  # Transcriptional misregulation in cancer
    'ko05203',  # Viral carcinogenesis
    'ko05020',  # Prion diseases
    'ko05143',  # African trypanosomiasis
    'ko05142',  # Chagas disease (American trypanosomiasis)
    'ko05140',  # Leishmaniasis
    'ko00072',  # Synthesis and degradation of ketone bodies
    'ko00073',  # Cutin, suberine and wax biosynthesis
    'ko05330',  # Allograft rejection
    'ko04912',  # GnRH signaling pathway
    'ko05332',  # Graft-versus-host disease
    'ko04910',  # Insulin signaling pathway
    'ko04916',  # Melanogenesis
    'ko04914',  # Progesterone-mediated oocyte maturation
    'ko04672',  # Intestinal immune network for IgA production
    'ko04670',  # Leukocyte transendothelial migration
    'ko04660',  # T cell receptor signaling pathway
    'ko00603',  # Glycosphingolipid biosynthesis - globo series
    'ko00601',  # Glycosphingolipid biosynthesis - lacto and neolacto series
    'ko00600',  # Sphingolipid metabolism
    'ko00604',  # Glycosphingolipid biosynthesis - ganglio series
    'ko04070',  # Phosphatidylinositol signaling system
    'ko04075',  # Plant hormone signal transduction
    'ko05100',  # Bacterial invasion of epithelial cells
    'ko05204',  # Chemical carcinogenesis
    'ko04920',  # Adipocytokine signaling pathway
    'ko04080',  # Neuroactive ligand-receptor interaction
    'ko04320',  # Dorso-ventral axis formation
    'ko04730',  # Long-term depression
    'ko04950',  # Maturity onset diabetes of the young
    'ko04630',  # Jak-STAT signaling pathway
    'ko04742',  # Taste transduction
    'ko04740',  # Olfactory transduction
    'ko04744',  # Phototransduction
    'ko04745',  # Phototransduction - fly
    'ko04512',  # ECM-receptor interaction
    'ko03460',  # Fanconi anemia pathway
    'ko04360',  # Axon guidance
    'ko05219',  # Bladder cancer
    'ko05218',  # Melanoma
    'ko05217',  # Basal cell carcinoma
    'ko05216',  # Thyroid cancer
    'ko05215',  # Prostate cancer
    'ko05214',  # Glioma
    'ko05213',  # Endometrial cancer
    'ko05212',  # Pancreatic cancer
    'ko05211',  # Renal cell carcinoma
    'ko05210',  # Colorectal cancer
    'ko05034',  # Alcoholism
    'ko05033',  # Nicotine addiction
    'ko05032',  # Morphine addiction
    'ko05031',  # Amphetamine addiction
    'ko05030',  # Cocaine addiction
    'ko05134',  # Legionellosis
    'ko05132',  # Salmonella infection
    'ko05133',  # Pertussis
    'ko05130',  # Pathogenic Escherichia coli infection
    'ko05131',  # Shigellosis

    #troppo grandi
    'ko01100',  # Metabolic pathways
    'ko01110',  # Biosynthesis of secondary metabolites
    'ko01120',  # Microbial metabolism in diverse environments

    'ko04662',  # B cell receptor signaling pathway
    'ko04966',  # Collecting duct acid secretion
    'ko05322',  # Systemic lupus erythematosus
    'ko04964',  # Proximal tubule bicarbonate reclamation
    'ko05320',  # Autoimmune thyroid disease
    'ko04962',  # Vasopressin-regulated water reabsorption
    'ko04960',  # Aldosterone-regulated sodium reabsorption
    'ko04961',  # Endocrine and other factor-regulated calcium reabsorption
    'ko04380',  # Osteoclast differentiation
    'ko04664',  # Fc epsilon RI signaling pathway
    'ko04666',  # Fc gamma R-mediated phagocytosis
    'ko04062',  # Chemokine signaling pathway
    'ko04060',  # Cytokine-cytokine receptor interaction
    'ko04066',  # HIF-1 signaling pathway
    'ko04064',  # NF-kappa B signaling pathway
    'ko05323',  # Rheumatoid arthritis
    'ko00140',  # Steroid hormone biosynthesis
]


def download_data(fname='kegg.pickle', contact=None):

    kclient = KeggClientRest()
    kclient.contact = contact
    kdata = KeggData()

    LOG.info("Downloading pathway list")
    path_names = kclient.get_ortholog_pathways()
    LOG.info("Found %d pathways", len(path_names))

    LOG.info("Downloading KO descriptions")
    ko_names = kclient.get_kos_descriptions()
    LOG.info("Found %d KOs", len(ko_names))

    LOG.info("Downloading reactions descriptions")
    rn_names = kclient.get_reactions_descriptions()
    LOG.info("Found %d reactions", len(rn_names))

    LOG.info("Downloading compounds descriptions")
    cpd_names = kclient.get_compounds_descriptions()
    LOG.info("Found %d compounds", len(cpd_names))

    kos = {}
    rns = {}
    cps = {}

    LOG.info("Downloading links pathway-ko (%d)", len(path_names))
    path_links = kclient.link_ids('ko', path_names.keys())
    for path_id, ko_list in path_links.iteritems():
        path = KeggPathway(path_id, path_names[path_id])
        for ko_id in ko_list:
            try:
                name = ko_names[ko_id]
            except KeyError:
                #in this case the actual entry for this gene doesn't exists
                LOG.warning(
                    "KO %s not found in the descriptions, skipping", ko_id
                )
                continue
            try:
                ko = kos[ko_id]
            except KeyError:
                ko = KeggOrtholog(ko_id, name)
                kos[ko_id] = ko
            path[ko_id] = ko
        kdata.pathways[path_id] = path

    LOG.info("Downloading links ko-reactions (%d)", len(ko_names))
    ko_links = kclient.link_ids('rn', ko_names.keys())
    for ko_id, rn_list in ko_links.iteritems():
        try:
            ko = kos[ko_id]
        except KeyError:
            LOG.warning("KO %s is not found in a pathway, skipping", ko_id)
            continue
        for rn_id in rn_list:
            try:
                name = rn_names[rn_id]
            except KeyError:
                LOG.warning(
                    "Reaction %s not found in the descriptions," +
                    " skipping", rn_id)
            try:
                rn = rns[rn_id]
            except KeyError:
                rn = KeggReaction(rn_id, name)
                rns[rn_id] = rn
            ko[rn_id] = rn

    LOG.info("Downloading links reactions-compounds (%d)", len(rn_names))
    cp_links = kclient.get_reaction_equations(rn_names.keys())
    for rn_id, cp_dict in cp_links.iteritems():
        try:
            rn = rns[rn_id]
        except KeyError:
            LOG.warning("Reaction %s is not found in a pathway, skipping",
                        rn_id)
            continue
        for cp_id in cp_dict['in'] + cp_dict['out']:
            try:
                name = cpd_names[cp_id]
            except KeyError:
                LOG.warning("Compound %s not found in the descriptions, " +
                            "skipping", rn_id)
            try:
                cp = cps[cp_id]
            except KeyError:
                cp = KeggCompound(cp_id, name)
                cps[cp_id] = cp
            if (cp_id in cp_dict['in']) and (cp_id in cp_dict['out']):
                LOG.debug(
                    "Compound %s in both side of reaction %s",
                    cp_id, rn_id
                )
            if cp_id in cp_dict['in']:
                rn.cp_in[cp_id] = cp
            if cp_id in cp_dict['out']:
                rn.cp_out[cp_id] = cp

    kdata.save_data(fname)


class KeggModule(object):
    """
    .. versionadded:: 0.1.13

    Used to extract information from a pathway module entry in Kegg

    The entry, as a string, can be either passed at instance creation or with
    :meth:`KeggModule.parse_entry`

    """
    entry = ''
    name = ''
    classes = None
    compounds = None

    _orthologs = None
    _reactions = None
    reactions = None

    def __init__(self, entry=None):
        if entry is None:
            return
        self.parse_entry(entry)

    def parse_entry(self, entry):
        """
        Parses a Kegg module entry and change the instance values. By default
        the reactions IDs are substituted with the KO IDs
        """
        entryd = {}
        curr_field = ''
        for line in entry.splitlines():
            if line.startswith(' '):
                entryd[curr_field].append(line.strip())
            elif line.startswith('///'):
                continue
            else:
                curr_field = line.split(' ')[0]
                entryd[curr_field] = []
                entryd[curr_field].append(line.replace(curr_field, '').strip())

        self.entry = re.search(r"(M\d{5})\s+.+", entryd['ENTRY'][0]).group(1)
        self.name = entryd['NAME'][0]

        self.classes = entryd['CLASS'][0].split('; ')
        self.compounds = [
            re.search(r"(C\d{5})\s+.+", line).group(1)
            for line in entryd['COMPOUND']
        ]

        self.reactions = [
            self.parse_reaction(reaction, ko_ids)
            for ko_ids, reaction in zip(entryd['DEFINITION'][0].split(' '), entryd['REACTION'])
        ]
        self._orthologs = entryd['DEFINITION'][0].split(' ')
        self._reactions = entryd['REACTION']

    @staticmethod
    def parse_reaction(line, ko_ids=None):
        """
        parses the lines with the reactions and substitute reaction IDs with the
        corresponding KO IDs if provided
        """

        #line = 'R00294,R02492,R09446,R09808,R09809  C00533 -> C00887'
        #ko_ids = '(K00370+K00371+K00374+K00373,K02567+K02568)'
        #some reaction lines have only one space
        rn_ids, reaction = line.replace('  ', ' ').split(' ', 1)
        rn_ids = rn_ids.split(',')
        comp1, comp2 = reaction.split(' -> ')
        if ko_ids is not None:
            rn_ids = ko_ids.replace('+', ',').replace('-', ',').replace('(', '').replace(')', '').split(',')
        return tuple(rn_ids), (comp1, comp2)

    @property
    def first_cp(self):
        "Returns the first compound in the module"
        return self.compounds[0]

    @property
    def last_cp(self):
        "Returns the first compound in the module"
        return self.compounds[-1]

    def to_edges(self):
        "Returns the reactions as edges that can be supplied to a graph"
        for ko_ids, (comp1, comp2) in self.reactions:
            for ko_id in ko_ids:
                yield (comp1, ko_id)
                yield (ko_id, comp2)
