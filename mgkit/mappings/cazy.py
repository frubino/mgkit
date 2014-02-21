"""
Module containing classes and functions to deal with CaZy data
"""

import urllib2
import logging
from .. import kegg

CAZY_FAMILIES = {
    'GH': 'Glycoside Hydrolase',
    'GT': 'GlycosylTransferase',
    'PL': 'Polysaccharide Lyase',
    'CE': 'Carbohydrate Esterase',
    'CBM': 'Carbohydrate-Binding Module'
}
"CaZy families"

LOG = logging.getLogger(__name__)


class Kegg2CazyMapper(kegg.KeggMapperBase):
    """
    Class holding mappings KOs->CaZy
    """
    query_string = 'database:(type:ko {0}) AND database:(type:cazy)'
    columns_string = 'database(cazy)'

    def map_kos_cazy(self, kos, contact, skip=False):
        """
        Map a list of KOs to CaZy ids
        """
        msg1 = "(%05d/%05d) Already downloaded mapping for KO: %s - %s"
        msg2 = "(%05d/%05d) Getting Cazy mappings for KO: %s - %s"
        warn1 = "(%05d/%05d) Poblem downloading KO %s (%s)"

        errors = 0
        for idx, ko_id in enumerate(kos):
            if skip and (ko_id in self._ko_map or ko_id in self._not_found):
                LOG.info(msg1, idx + 1, len(kos), ko_id, kos[ko_id])
                continue
            LOG.info(msg2, idx + 1, len(kos), ko_id, kos[ko_id])
            try:
                mappings = self.ko_to_mapping(
                    ko_id, self.query_string,
                    self.columns_string, contact
                )
            except urllib2.HTTPError as url_error:
                LOG.warning(warn1, idx + 1, len(kos), ko_id, url_error.args)
                errors += 1
                continue
            if mappings:
                self._ko_map[ko_id] = mappings
            else:
                self._not_found.append(ko_id)
        return errors


def download_data(contact, kegg_data='kegg.pickle', cazy_data='cazy.pickle'):
    """
    Function to download CaZy data
    """
    WARN1 = "Cazy data not found, download restarting"
    MSG1 = "Loading Kegg data from file %s"
    MSG2 = "Saving downloaded data to %s"
    MSG3 = "Found %d/%d mappings. Number of errors: %d"

    c_obj = Kegg2CazyMapper()
    try:
        c_obj.load_data(cazy_data)
    except IOError:
        LOG.warning(WARN1)

    LOG.info(MSG1, kegg_data)
    kegg_obj = kegg.KeggData(kegg_data)
    ko_names = kegg_obj.get_ko_names()

    try:
        errors = c_obj.map_kos_cazy(ko_names, contact, skip=True)
    finally:
        LOG.info(MSG2, cazy_data)
        c_obj.save_data(cazy_data)

    LOG.info(MSG3, len(c_obj), len(ko_names), errors)
