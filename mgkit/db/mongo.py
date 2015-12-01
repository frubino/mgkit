"""
.. versionadded:: 0.2.1

This module contains functions and classes to use for a DB like representation
of annotations using the *pymongo* package, a driver to **MongoDB**.

In a MongoDB document, exported from an annotation, using the
:meth:`mgkit.io.gff.Annotation.to_mongodb` method, the keys that are defined
are::

    seq_id, source, feat_type, start, end, score, strand,
    phase, gene_id, taxon_id, bitscore, exp_nonsyn, exp_syn,
    length, dbq, coverage, map

These are defined because they have values that are not strings (defined as
properties in :class:`mgkit.io.gff.Annotation`. The rest of the attributes
defined are kept as well, but no ckeck for the data type is made.

.. note::

    lineage is added as a key, whose values are taxon_id, if a function has been
    passed to :meth:`mgkit.io.gff.Annotation.to_mongodb`

The exception is the **map** key in the document. It store both the EC mappings
(EC attribute in the GFF), as well as all mappings whose attribute starts with
*map_*. The former is usually accessed from
:meth:`mgkit.io.gff.Annotation.get_ec` while the latter from
:meth:`mgkit.io.gff.Annotation.get_mapping` or
:meth:`mgkit.io.gff.Annotation.get_mappings`.

These 3 methods return a list and this list is used in the MongoDB document.
The MongoDB document will contain a **map** key where the values are the type
of mappings, and the values the list of IDs the annoation maps to.

.. list-table:: Example for the map dictionary
   :header-rows: 1

   * - Type
     - GFF
     - Annotation
     - MongoDB Document
     - MongoDB Query
   * - EC
     - EC
     - get_ec
     - ec
     - map.ec
   * - KO
     - map_KO
     - get_mapping('ko')
     - ko
     - map.ko
   * - eggNOG
     - map_EGGNOG
     - get_mapping('eggnog')
     - eggnog
     - map.eggnog


"""
import logging
from pymongo import MongoClient
from ..io import gff

LOG = logging.getLogger(__name__)


class GFFDB(object):
    """
    Wrapper to a MongoDB connection/db. It is used to automate the convertion
    of MongoDB records into :class:`mgkit.io.gff.Annotation` instances.
    """
    conn = None
    db = None

    def __init__(self, db, collection, uri=None):
        self.conn = MongoClient(uri)
        self.db = self.conn[db][collection]

    def cursor(self, query=None):
        "Returns a cursor for the query"
        return self.db.find(query)

    def convert_record(self, record):
        "Converts the record (a dictionary instance) to an Annotation"
        return gff.from_mongodb(record)

    def find_annotation(self, query=None):
        """
        Iterate over a cursor created using *query* and yields each record
        after converting it to a :class:`mgkit.io.gff.Annotation` instance,
        using :meth:`mgkit.db.mongo.GFFDB.convert_record`.
        """
        for record in self.cursor(query):
            yield self.convert_record(record)
