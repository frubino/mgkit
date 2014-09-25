"""
Package used to contain code related to I/O operations
"""

import gzip
import bz2

try:
    import lzma
except ImportError:
    lzma = None


class UnsupportedFormat(IOError):
    "Raised if the a file can't be opened with the correct module"
    pass


def open_file(file_name, mode='r'):
    """
    .. versionadded:: 0.1.12

    Opens a file using the extension as a guide to which module to use.

    Arguments:
        file_name (str): file name
        mode (str): mode used to open the file

    Returns:
        file: file handle

    Raises:
        UnsupportedFormat: if the module to open the file is not available

    """
    if file_name.endswith('.gz'):
        file_handle = gzip.GzipFile(file_name, mode)
    elif file_name.endswith('.bz2'):
        file_handle = bz2.BZ2File(file_name, mode)
    elif file_name.endswith('.xz'):
        if lzma:
            raise UnsupportedFormat("Cannot import lzma module")
        else:
            file_handle = lzma.LZMAFile(file_name, mode)
    else:
        file_handle = open(file_name, mode)

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
    out_handles = [open(name_mask.format(x), 'w') for x in xrange(num_files)]

    for index, record in enumerate(records):
        out_handle = out_handles[index % num_files]
        write_func(out_handle, record)
