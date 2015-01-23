# coding=utf8
"""
Simple fasta parser and a few utility functions
"""
import logging
import mgkit.io

LOG = logging.getLogger(__name__)


def load_fasta(f_handle):
    """
    .. versionchanged:: 0.1.13
        now returns uppercase sequences

    Loads a fasta file and returns a generator of tuples in which the first
    element is the name of the sequence and the second the sequence

    Arguments:
        f_handle (str, file): fasta file to open; a file name or a file handle
            is expected

    Yields:
        tuple: first element is the sequence name/header, the second element is
        the sequence
    """
    if isinstance(f_handle, str):
        f_handle = mgkit.io.open_file(f_handle, 'r')

    if getattr(f_handle, 'name', None) is not None:
        LOG.info("Reading fasta file %s", f_handle.name)

    cur_name = ""
    last_name = ""
    nseq = 0
    # Better use ''.join(list) than seq+=seq, it's faster
    cur_seq = []
    # main loop to read file's sequences
    for line in f_handle:
        if line.startswith('>'):
            if cur_seq != []:
                # start of next sequence
                yield cur_name, ''.join(cur_seq).upper()
                nseq += 1
                cur_seq = []
            # save previous name for loop's else clause
            last_name = cur_name
            # classic fasta, all characters besides starting ">" and rightmost
            # whitespace
            cur_name = line[1:].rstrip()
        else:
            cur_seq.append(line.rstrip())
    else:
        # handles cases where there's no newline after last
        # portion of the sequence at EOF
        if nseq != 0:
            if cur_name != last_name:
                yield cur_name, ''.join(cur_seq).upper()
        # case in which only one sequence is present
        else:
            yield cur_name, ''.join(cur_seq).upper()
    f_handle.close()


def write_fasta_sequence(file_handle, name, seq, wrap=60, write_mode='a'):
    """
    Write a fasta sequence to file. If the *file_handle* is a string, the file
    will be opened using *write_mode*.

    :param file_handle: file handle or string.
    :param str name: header to write for the sequence
    :param str seq: sequence to write
    :param int wrap: int for the line wrapping. If None, the sequence will be
        written in a single line
    """
    if isinstance(file_handle, str):
        file_handle = open(file_handle, write_mode)

    if wrap is not None:
        seq = '\n'.join(
            seq[pos:pos + wrap]
            for pos in xrange(0, len(seq), wrap)
        )

    file_handle.write(">{0}\n{1}\n".format(name, seq))


def split_fasta_file(file_handle, name_mask, num_files):
    """
    .. versionadded:: 0.1.13

    Splits a fasta file into a series of smaller files.

    Arguments:
        file_handle (file, str): fasta file with the input sequences
        name_mask (str): file name template for the splitted files, more
            informations are found in :func:`mgkit.io.split_write`
        num_files (int): number of files in which to distribute the sequences
    """
    records = load_fasta(file_handle)

    write_func = lambda h, r: write_fasta_sequence(h, *r)

    mgkit.io.split_write(records, name_mask, write_func, num_files=num_files)
