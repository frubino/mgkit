"""
Various utilities to help read and process files
"""
from builtins import range, next
import sys
import logging
import gzip
import bz2
import io

try:
    if sys.version_info >= (3, 3):
        import lzma
    else:
        from backports import lzma
except ImportError:
    lzma = None

LOG = logging.getLogger(__name__)


def group_tuples_by_key(iterator, key_func=None, skip_elements=0):
    """
    .. versionadded:: 0.3.1

    Group the elements of an iterator by a key and yields the grouped elements.
    The elements yielded by the iterator are assumed to be a list or tuple,
    with the default key (when *key_func* is None) being the first of the of
    the objects inside that element. This behaviour can be customised by
    passing to *key_func* a function that accept an element and returns the key
    to be used.

    .. note::

        the iterable assumen that the elements are already sorted by their keys

    Arguments:
        iterator (iterable): iterator to be grouped
        key_func (func): function that accepts a element and returns its
            associated key
        skip_elements (int): number of elements to skip at the start

    Yields:
        list: a list of the grouped elements by key
    """
    if key_func is None:
        def key_func(x): return x[0]

    for index in range(skip_elements):
        next(iterator)

    curr_key = None
    curr_ann = []

    for element in iterator:
        new_key = key_func(element)
        if curr_key == new_key:
            curr_ann.append(element)
        else:
            if curr_key is None:
                curr_ann.append(element)
                curr_key = new_key
            else:
                yield curr_ann
                curr_key = new_key
                curr_ann = [element]
    else:
        yield curr_ann


class UnsupportedFormat(IOError):
    "Raised if the a file can't be opened with the correct module"
    pass


def open_file(file_name, mode='r'):
    """
    .. versionadded:: 0.1.12

    .. versionchanged:: 0.3.4
        using *io.open*, always in binary mode

    .. versionchanged:: 0.4.2
        when a file handle is detected, it is passed to :func:`compressed_handle`
        to detect if the handle is a compressed file

    Opens a file using the extension as a guide to which module to use.

    .. note::

        Unicode makes for a slower `.translate` method in Python2, so it's
        best to use the `open` builtin.

    Arguments:
        file_name (str): file name
        mode (str): mode used to open the file

    Returns:
        file: file handle

    Raises:
        UnsupportedFormat: if the module to open the file is not available

    """

    if sys.version_info[0] == 2:
        test_class = (file, io.IOBase)
    else:
        test_class = io.IOBase

    if isinstance(file_name, test_class):
        return compressed_handle(file_name)

    mode = mode + 'b' if 'b' not in mode else mode

    if file_name.endswith('.gz'):
        file_handle = gzip.GzipFile(file_name, mode)
    elif file_name.endswith('.bz2'):
        file_handle = bz2.BZ2File(file_name, mode)
    elif file_name.endswith('.xz'):
        if lzma is None:
            raise UnsupportedFormat("Cannot import lzma module")
        else:
            file_handle = lzma.LZMAFile(file_name, mode)
    else:
        file_handle = io.open(file_name, mode)

    return file_handle


def compressed_handle(file_handle):
    """
    .. versionadded:: 0.1.13

    Tries to wrap a file handle in the appropriate compressed file class.

    Arguments:
        file_handle (str): file handle

    Returns:
        file: the same file handle if no suitable compressed file class is
        found or the new file_handle which supports the compression

    Raises:
        UnsupportedFormat: if the module to open the file is not available

    """
    LOG.info('Detecting Compressed File Handle')

    if file_handle.name.endswith('.gz'):
        file_handle = gzip.GzipFile(fileobj=file_handle, mode='rb')
    elif file_handle.name.endswith('.xz'):
        if lzma:
            raise UnsupportedFormat("Cannot import lzma module")
        else:
            file_handle = lzma.LZMAFile(file_handle)

    return file_handle


def split_write(records, name_mask, write_func, num_files=2):
    """
    .. versionadded:: 0.1.13

    Splits the writing of a number of records in a series of files. The
    `name_mask` is used as template for the file names. A string like
    "split-files-{0}" can be specified and the function applies format with the
    index of the pieces.

    Arguments:
        records (iterable): an iterable that returns a object to be saved
        name_mask (str): a string used as template for the output file names
            on which the function applies :func:`string.format`
        write_func (func): a function that is called to write to the files. It
            needs to accept a file handles as first argument and the record
            returned by `records` as the second argument
        num_files (int): the number of files to split the records
    """
    out_handles = [open_file(name_mask.format(x), 'w') for x in range(num_files)]

    for index, record in enumerate(records):
        out_handle = out_handles[index % num_files]
        write_func(out_handle, record)
