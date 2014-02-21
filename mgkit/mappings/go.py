
"""
Module containing classes and functions to deal with Gene Ontology data
"""

import urllib2
import logging
from .. import kegg
from ..utils import dictionary as dict_utils
from goatools.obo_parser import OBOReader

LOG = logging.getLogger(__name__)


class Kegg2GOMapper(kegg.KeggMapperBase):
    """
    Class holding mappings KOs->GO
    """
    query_string = 'database:(type:ko {0})'
    columns_string = 'go-id'
    _id_names = {}

    def map_kos_go(self, kos, contact, skip=False):
        """
        Map a list of KOs to GO ids
        """
        msg1 = "(%05d/%05d) Already downloaded mapping for KO: %s - %s"
        msg2 = "(%05d/%05d) Getting GO mappings for KO: %s - %s"
        msg3 = "(%05d/%05d) Poblem downloading KO %s (%s)"

        errors = 0
        for idx, ko_id in enumerate(kos):
            if skip and (ko_id in self._ko_map or ko_id in self._not_found):
                LOG.info(msg1, idx + 1, len(kos), ko_id, kos[ko_id])
                continue
            LOG.info(msg2, idx + 1, len(kos), ko_id, kos[ko_id])
            try:
                mappings = self.ko_to_mapping(ko_id, self.query_string,
                                              self.columns_string, contact)
            except urllib2.HTTPError as url_error:
                LOG.warning(msg3, idx + 1, len(kos), ko_id, url_error.args)
                errors += 1
                continue
            if mappings:
                self._ko_map[ko_id] = mappings
            else:
                self._not_found.append(ko_id)
        return errors

    def load_goslim_mapping(self, f_handle):
        if isinstance(f_handle, str):
            f_handle = open(f_handle, 'r')

        goslim_map = {}

        for line in f_handle:
            line = line.split('//')[0]
            goterm, goslim = line.split('=>')
            goterm = goterm.strip()
            if not goterm.startswith('GO'):
                continue
            goslim = goslim.strip().split()
            goslim_map[goterm] = goslim

        self._ko_map = dict_utils.combine_dict(self._ko_map, goslim_map)

    def load_go_names(self, fname):
        name_map = self._id_names
        for term in OBOReader(fname):
            name_map[term.id] = term.name

    def write_association_file(self, f_handle):
        if isinstance(f_handle, str):
            f_handle = open(f_handle, 'w')

        ko_map = self.get_ko_map()

        for ko_id, goterms in ko_map.iteritems():
            f_handle.write("{0}\t{1}\n".format(ko_id, ';'.join(goterms)))


def download_data(contact, kegg_data='kegg.pickle', go_data='go_data.pickle'):
    """
    Function to download Gene Ontology data
    """

    g_obj = Kegg2GOMapper()
    try:
        g_obj.load_data(go_data)
    except IOError:
        LOG.warning("GO data not found, download restarting")

    LOG.info("Loading Kegg data from file %s", kegg_data)
    kegg_obj = kegg.KeggData(kegg_data)
    ko_names = kegg_obj.get_ko_names()

    try:
        errors = g_obj.map_kos_go(ko_names, contact, skip=True)
    finally:
        LOG.info("Saving downloaded data to %s", go_data)
        g_obj.save_data(go_data)

    LOG.info(
        "Found %d/%d mappings. Number of errors: %d",
        len(g_obj), len(ko_names), errors
    )
