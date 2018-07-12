"""
.. versionadded:: 0.2.1

This module contains functions and classes to use for a dbm like representation
of annotations using the *semidbm* package

"""
import logging
from builtins import object, bytes
from ..io import gff

import semidbm

LOG = logging.getLogger(__name__)


def create_gff_dbm(annotations, file_name):
    """
    .. versionadded:: 0.2.1

    Creates a semidbm database, using an annotation `uid` as key and the gff
    line as value. The object is synced before being returned.

    .. note::

        A GFF line is used instead of a json representation because it was
        more compact when semidbm was tested.

    Arguments:
        annotations (iterable): iterable of annotations
        file_name (str): database file name, opened with the `c` flag.

    Returns:
        object: a semidbm database object
    """
    database = semidbm.open(file_name, 'c')

    LOG.info('DB "%s" opened/created', file_name)

    for annotation in annotations:
        database[annotation.uid.encode('ascii')] = \
            annotation.to_gff().encode('ascii')

    database.sync()

    return database


class GFFDB(object):
    """
    .. versionadded:: 0.2.1

    A wrapper for a semidbm instance, used to convert the GFF line stored in
    the DB into an :class:`mgkit.io.gff.Annotation` instance. If a string is
    passed to the init method, a DB will be opened with the `c` flag.

    The object behaves like a dictionary, wrapping the access to annoations
    using a *uid* as key and converting the line into an
    :class:`mgkit.io.gff.Annotation` instance.
    """
    db = None

    def __init__(self, db=None):
        if isinstance(db, str):
            self.db = semidbm.open(db, 'c')
        else:
            self.db = db

    def __setitem__(self, uid, annotation):
        self.db[uid.encode('ascii')] = annotation.to_gff().encode('ascii')

    def __getitem__(self, key):
        if not isinstance(key, bytes):
            key = key.encode('ascii')
        return gff.from_gff(self.db[key].decode('ascii'))

    def __del__(self):
        self.db.close()

    def __iter__(self):
        for uid in self.db:
            yield uid.decode('ascii')

    def items(self):
        for uid in self:
            yield uid, self[uid]

    def iteritems(self):
        return self.items()

    def values(self):
        for uid in self:
            yield self[uid]

    def itervalues(self):
        return self.values()
