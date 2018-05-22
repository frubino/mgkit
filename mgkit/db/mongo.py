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

    lineage is added as a key, whose values are taxon_id, if a function has
    been passed to :meth:`mgkit.io.gff.Annotation.to_mongodb`

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
from builtins import object
import logging
from ..io import gff

from pymongo import MongoClient

LOG = logging.getLogger(__name__)


class GFFDB(object):
    """
    .. versionchanged:: 0.3.4
        added *timeout* parameter

    Wrapper to a MongoDB connection/db. It is used to automate the convertion
    of MongoDB records into :class:`mgkit.io.gff.Annotation` instances.
    """
    conn = None
    db = None

    def __init__(self, db, collection, uri=None, timeout=5):
        self.conn = MongoClient(uri, serverSelectionTimeoutMS=timeout)
        self.db = self.conn[db][collection]

    def cursor(self, query=None):
        "Returns a cursor for the query"
        return self.db.find(query)

    def convert_record(self, record):
        """
        .. versionchanged:: 0.3.1
            removes *lineage* from the attributes

        Converts the record (a dictionary instance) to an Annotation
        """
        return gff.from_mongodb(record, lineage=False)

    def find_annotation(self, query=None):
        """
        Iterate over a cursor created using *query* and yields each record
        after converting it to a :class:`mgkit.io.gff.Annotation` instance,
        using :meth:`mgkit.db.mongo.GFFDB.convert_record`.
        """
        for record in self.cursor(query):
            yield self.convert_record(record)

    def __getitem__(self, uid):
        """
        .. versionadded:: 0.3.1

        Retrieves an annotation from the DB by its *uid*
        """
        return self.convert_record(self.db.find_one(uid))

    def insert_one(self, annotation):
        """
        .. versionadded:: 0.3.4

        Inserts an annotation into the DB

        Raises:
            TypeError: if the passed object is not an annotation
        """
        if not isinstance(annotation, gff.Annotation):
            raise TypeError(
                'An mgkit.io.gff.Annotation instance is required got {}'.format(
                    type(annotation)
                )
            )

        self.db.insert_one(annotation.to_mongodb(raw=True))

    def insert_many(self, annotations):
        """
        .. versionadded:: 0.3.4

        Inserts annotations into the DB

        .. warning::

            The object must be a :class:`mgkit.io.gff.Annotation`

        """

        self.db.insert_many(
            (
                annotation.to_mongodb(raw=True)
                for annotation in annotations
            ),
            ordered=False
        )

    def __iter__(self):
        """
        .. versionadded:: 0.3.1

        Iterates over all annotations
        """
        return self.values()

    def values(self):
        """
        .. versionadded:: 0.3.1

        Iterates over all the annotations in the db/collection
        """
        return self.find_annotation()

    def itervalues(self):
        """
        .. versionadded:: 0.3.1

        Alias for :meth:`GFFDB.values`
        """
        return self.values()

    def items(self):
        """
        .. versionadded:: 0.3.1

        Iterates over all the annotations in the db/collection, yielding a
        tuple (*annotation.uid*, *annotation*)
        """
        for record in self.values():
            yield record.uid, record

    def iteritems(self):
        """
        .. versionadded:: 0.3.1

        Alias for :meth:`GFFDB.items`
        """
        return self.items()

    def keys(self):
        """
        .. versionadded:: 0.3.1

        Iterates over all the *uid* in the db/collection
        """
        for record in self.values():
            yield record.uid
