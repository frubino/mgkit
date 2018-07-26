"""
Module containing classes and functions to access Kegg data
"""
from builtins import object
from future.utils import viewitems
import sys
import logging
import pickle
import random
import re
from .utils import dictionary as dict_utils
from .net import url_read
from .io import open_file


LOG = logging.getLogger(__name__)
KEGG_REST_URL = 'http://rest.kegg.jp/'


class KeggClientRest(object):
    """
    .. versionchanged:: 0.3.1
        added a *cache* attribute for some methods

    Kegg REST client

    The class includes methods and data to use the REST API provided by Kegg.
    At the moment it provides methods to for 'link', 'list' and 'get'
    operations,

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

    cache = None

    def __init__(self, cache=None):
        """
        .. versionadded:: 0.3.1

        The "cache" parameter is a file name for the cached data wrote using
        :meth:`KeggClientRest.write_cache`.
        """

        if cache is None:
            self.empty_cache()
        else:
            self.load_cache(cache)

    def empty_cache(self, methods=None):
        """
        .. versionadded:: 0.3.1

        Empties the cache completely or for a specific method(s)

        Arguments:
            methods (iterable, str): string or iterable of strings that are
                part of the cache. If None the cache is fully emptied
        """

        if methods is None:
            methods = ('link_ids', 'get_entry', 'get_ids_names')
        else:
            if isinstance(methods, str):
                methods = [methods]

        if self.cache is None:
            self.cache = {}

        for method in methods:
            self.cache[method] = {}

    def load_cache(self, file_handle):
        """
        .. versionadded:: 0.3.1

        Loads the cache from file
        """
        self.cache = pickle.load(open_file(file_handle, 'rb'))

    def write_cache(self, file_handle):
        """
        .. versionadded:: 0.3.1

        Write the cache to file
        """
        pickle.dump(self.cache, open_file(file_handle, 'wb'))

    # Kegg primitives #

    def link_ids(self, target, kegg_ids, max_len=50):
        """
        .. versionchanged:: 0.3.1
            removed *strip* and cached the results

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
        if (sys.version_info[0] == 2) and isinstance(kegg_ids, unicode):
            kegg_ids = [kegg_ids]
        if isinstance(kegg_ids, str):
            kegg_ids = [kegg_ids]

        if isinstance(kegg_ids, (list, dict)):
            kegg_ids = set(kegg_ids)

        try:
            mapping = self.cache['link_ids'][target]
            LOG.debug('Cached Call')
        except KeyError:
            LOG.debug('Empty Cache')
            mapping = {}
            self.cache['link_ids'][target] = mapping

        download_ids = kegg_ids - set(mapping)
        LOG.debug(
            'Number of cached items (%d/%d)',
            len(kegg_ids) - len(download_ids),
            len(kegg_ids)
        )

        while download_ids:
            try:
                sample_ids = set(random.sample(download_ids, max_len))
            except ValueError:
                # download_ids size is less than max_len
                sample_ids = download_ids.copy()
            if len(kegg_ids) > max_len:
                LOG.info(
                    "Downloading links - %d out of %d",
                    len(kegg_ids) - len(download_ids), len(kegg_ids)
                )
            url = "{0}link/{1}/{2}".format(
                self.api_url, target, '+'.join(sample_ids)
            )
            data = url_read(url, agent=self.contact)
            data = data.split('\n')

            for line in data:
                if not line:
                    continue
                source, linked = line.split('\t')
                source = source.split(':')[1]
                linked = linked.split(':')[1]

                try:
                    mapping[source].append(linked)
                except KeyError:
                    mapping[source] = [linked]

            download_ids = download_ids - sample_ids

        # Empty (Not Found) elements, set to None
        for kegg_id in kegg_ids - set(mapping):
            mapping[kegg_id] = None

        return {
            kegg_id: list(value)
            for kegg_id, value in viewitems(mapping)
            if (kegg_id in kegg_ids) and (value is not None)
        }

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
        # leave out the last \n
        return data[:-1]

    def get_entry(self, k_id, option=None):
        """
        .. versionchanged:: 0.3.1
            this is now cached

        The method abstract the use of the 'get' operation in the Kegg API

        :param str k_id: kegg id of the resource to get
        :param str option: optional, to specify a format
        """
        try:
            data = self.cache['get_entry'][(k_id, option)]
        except KeyError:
            url = "{0}get/{1}/{2}".format(
                self.api_url,
                k_id,
                '' if option is None else option
            )
            data = url_read(url, agent=self.contact)
            self.cache['get_entry'][(k_id, option)] = data

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

    def find(self, query, database, options=None, strip=True):
        """
        .. versionadded:: 0.3.1

        Kegg Help:

        http://rest.kegg.jp/find/<database>/<query>

        <database> = pathway | module | ko | genome | <org> | compound | glycan |
                     reaction | rclass | enzyme | disease | drug | dgroup | environ |
                     genes | ligand

        <org> = KEGG organism code or T number

        http://rest.kegg.jp/find/<database>/<query>/<option>

        <database> = compound | drug
        <option> = formula | exact_mass | mol_weight

        Examples:
            >>> kc = KeggClientRest()
            >>> kc.find('CH4', 'compound')
            {'C01438': 'Methane; CH4'}
            >>> kc.find('K00844', 'genes', strip=False)
            {'tped:TPE_0072': 'hexokinase; K00844 hexokinase [EC:2.7.1.1]',
            ...
            >>> kc.find('174.05', 'compound', options='exact_mass')
            {'C00493': '174.052823',
             'C04236': '174.052823',
             'C16588': '174.052823',
             'C17696': '174.052823',
             'C18307': '174.052823',
             'C18312': '174.052823',
             'C21281': '174.052823'}
        """

        url = 'http://rest.kegg.jp/find/{}/{}/{}'.format(
            database,
            query,
            '' if options is None else options
        )

        LOG.debug(url)

        data = url_read(url, agent=self.contact)

        mappings = {}
        for line in data.splitlines():
            target_id, description = line.rstrip().split('\t')

            if strip:
                target_id = target_id.split(':')[1]

            mappings[target_id] = description

        return mappings

    def conv(self, target_db, source_db, strip=True):
        """
        .. versionadded:: 0.3.1

        Kegg Help:

        http://rest.kegg.jp/conv/<target_db>/<source_db>

        (<target_db> <source_db>) = (<kegg_db> <outside_db>) | (<outside_db> <kegg_db>)

        For gene identifiers:
        <kegg_db> = <org>
        <org> = KEGG organism code or T number
        <outside_db> = ncbi-proteinid | ncbi-geneid | uniprot

        For chemical substance identifiers:
        <kegg_db> = drug | compound | glycan
        <outside_db> = pubchem | chebi
        http://rest.kegg.jp/conv/<target_db>/<dbentries>

        For gene identifiers:
        <dbentries> = database entries involving the following <database>
        <database> = <org> | genes | ncbi-proteinid | ncbi-geneid | uniprot
        <org> = KEGG organism code or T number

        For chemical substance identifiers:
        <dbentries> = database entries involving the following <database>
        <database> = drug | compound | glycan | pubchem | chebi

        Examples:
            >>> kc = KeggClientRest()
            >>> kc.conv('ncbi-geneid', 'eco')
            {'eco:b0217': {'ncbi-geneid:949009'},
             'eco:b0216': {'ncbi-geneid:947541'},
             'eco:b0215': {'ncbi-geneid:946441'},
             'eco:b0214': {'ncbi-geneid:946955'},
             'eco:b0213': {'ncbi-geneid:944903'},
            ...
            >>> kc.conv('ncbi-proteinid', 'hsa:10458+ece:Z5100')
            {'10458': {'NP_059345'}, 'Z5100': {'AAG58814'}}
        """

        url = 'http://rest.kegg.jp/conv/{}/{}/'.format(
            target_db,
            source_db,
        )

        LOG.debug(url)

        data = url_read(url, agent=self.contact)

        mappings = {}
        for line in data.splitlines():
            source_id, target_id = line.rstrip().split('\t')

            if strip:
                target_id = target_id.split(':')[1]
                source_id = source_id.split(':')[1]

            try:
                mappings[source_id].add(target_id)
            except KeyError:
                mappings[source_id] = set([target_id])

        return mappings

    # names #

    def get_ids_names(self, target='ko', strip=True):
        """
        .. versionadded:: 0.1.13

        .. versionchanged:: 0.3.1
            the call is now cached

        Returns a dictionary with the names/description of all the id of a
        specific target, (ko, path, cpd, etc.)

        If strip=True the id will stripped of the module abbreviation (e.g.
        md:M00002->M00002)
        """

        if strip:
            try:
                return self.cache['get_ids_names'][target].copy()
            except KeyError:
                LOG.debug('No cached values for "%s"', target)

        id_names = {}

        for line in self.list_ids(target).splitlines():
            kegg_id, name = line.strip().split('\t')
            if strip:
                kegg_id = kegg_id.split(':')[1]
            id_names[kegg_id] = name

        if strip:
            self.cache['get_ids_names'][target] = id_names.copy()

        return id_names

    def get_ortholog_pathways(self):
        """
        Gets ortholog pathways, replace 'map' with 'ko' in the id
        """
        data = self.get_ids_names('pathway')
        pathways = {}
        for kegg_id, name in viewitems(data):

            kegg_id = kegg_id.replace('map', 'ko')
            pathways[kegg_id] = name

        return pathways

    # end names #

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

    # Too big
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

    def __init__(self, entry=None, old=False):
        """
        .. versionchanged:: 0.3.0
            added *old* parameter, to use the old parser
        """
        if entry is None:
            return
        if old:
            self.parse_entry(entry)
        else:
            self.parse_entry2(entry)

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
            for ko_ids, reaction in zip(
                entryd['DEFINITION'][0].split(' '),
                entryd['REACTION']
            )
        ]
        self._orthologs = entryd['DEFINITION'][0].split(' ')
        self._reactions = entryd['REACTION']

    def parse_entry2(self, entry):
        """
        .. versionadded:: 0.3.0

        Parses a Kegg module entry and change the instance values. By default
        the reactions IDs are NOT substituted with the KO IDs.
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

        self.reactions = []

        for reaction in entryd['REACTION']:
            reaction = self.parse_reaction(reaction, ko_ids=None)
            self.reactions.append(reaction)

        self._orthologs = entryd['DEFINITION'][0].split(' ')
        self._reactions = entryd['REACTION']

    @staticmethod
    def parse_reaction(line, ko_ids=None):
        """
        .. versionchanged:: 0.3.0
            cleaned the parsing

        parses the lines with the reactions and substitute reaction IDs with
        the corresponding KO IDs if provided
        """

        # line = 'R00294,R02492,R09446,R09808,R09809  C00533 -> C00887'
        # ko_ids = '(K00370+K00371+K00374+K00373,K02567+K02568)'
        # some reaction lines have only one space
        rn_ids, reaction = line.replace('  ', ' ').split(' ', 1)
        rn_ids = tuple(x.strip() for x in rn_ids.replace('+', ',').split(','))
        comp1, comp2 = reaction.split(' -> ')
        comp1 = comp1.replace('(spontaneous)', '')
        comp2 = comp2.replace('(spontaneous)', '')
        comp1 = tuple(x.strip() for x in comp1.replace('+', ',').split(','))
        comp2 = tuple(x.strip() for x in comp2.replace('+', ',').split(','))
        if ko_ids is not None:
            rn_ids = ko_ids.replace('+', ',').replace('-', ',').replace('(', '').replace(')', '').split(',')
        return rn_ids, (comp1, comp2)

    @property
    def first_cp(self):
        "Returns the first compound in the module"
        return self.reactions[0][1][0][0]

    @property
    def last_cp(self):
        "Returns the first compound in the module"
        return self.reactions[-1][-1][-1][0]

    def to_edges(self, id_only=None):
        """
        .. versionchanged:: 0.3.0
            added id_only and changed to reflect changes in :attr:`reactions`

        Returns the reactions as edges that can be supplied to make graph.

        Arguments:
            id_only (None, iterable): if None the returned edges are for the
                whole module, if an iterable (converted to a :class:`set`),
                only edges for those reactions are returned

        Yield:
            tuple: the elements are the compounds and reactions in the module
        """

        if id_only is not None:
            id_only = set(id_only)

        for rn_ids, (comp1s, comp2s) in self.reactions:
            for rn_id in rn_ids:
                if (id_only is not None) and (rn_id not in id_only):
                    continue
                for comp1 in comp1s:
                    yield (comp1, rn_id)
                for comp2 in comp2s:
                    yield (rn_id, comp2)

    def find_submodules(self):
        """
        .. versionadded:: 0.3.0

        Returns the possible submodules, as a list of tuples where the elements
        are the first and last compounds in a submodule
        """
        sub_modules = []
        sub_module = None
        for rn_ids, (left_cpds, right_cpds) in self.reactions:
            if sub_module is None:
                sub_module = [left_cpds, right_cpds]
                continue

            if set(sub_module[1]) & set(left_cpds):
                sub_module[1] = right_cpds
            else:
                sub_modules.append((sub_module[0][0], sub_module[-1][-1]))
                sub_module = [left_cpds, right_cpds]
        else:
            sub_modules.append((sub_module[0][0], sub_module[-1][-1]))
        return sub_modules


def parse_reaction(line, prefix=('C', 'G')):
    """
    .. versionadded:: 0.3.1

    Parses a reaction equation from Kegg, returning the left and right
    components. Needs testing

    Arguments:
        line (str): reaction string

    Returns:
        tuple: left and right components as `sets`

    Raises:
        ValueError: if the

    """
    line = line.replace('EQUATION', '').strip()
    if '<=>' in line:
        line = line.replace(' ', '').split('<=>')
        left = set(x if x.startswith('C') else x[1:] for x in line[0].split('+'))
        right = set(x if x.startswith('C') else x[1:] for x in line[1].split('+'))
        return left, right
    elif '=>' in line:
        raise ValueError('>>>')
    elif '<=' in line:
        raise ValueError('<<<')
    else:
        raise ValueError('???')
