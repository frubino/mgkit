# coding=utf8
"""
Simple fasta parser and a few utility functions
"""
import itertools
import logging
from builtins import range
import mgkit.io

LOG = logging.getLogger(__name__)


def load_fasta(file_handle):
    """
    .. versionchanged:: 0.1.13
        now returns uppercase sequences

    Loads a fasta file and returns a generator of tuples in which the first
    element is the name of the sequence and the second the sequence

    Arguments:
        file_handle (str, file): fasta file to open; a file name or a file handle
            is expected

    Yields:
        tuple: first element is the sequence name/header, the second element is
        the sequence
    """
    file_handle = mgkit.io.open_file(file_handle, 'rb')

    if getattr(file_handle, 'name', None) is not None:
        LOG.info("Reading fasta file %s", file_handle.name)

    cur_name = ""
    last_name = ""
    nseq = 0
    # Better use ''.join(list) than seq+=seq, it's faster
    cur_seq = []
    # main loop to read file's sequences
    for line in file_handle:
        line = line.decode('ascii')
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
                nseq += 1
        # case in which only one sequence is present
        else:
            yield cur_name, ''.join(cur_seq).upper()
            nseq += 1
    file_handle.close()

    LOG.info("Read %d fasta sequences", nseq)


def load_fasta_files(files):
    """
    .. versionadded:: 0.3.4

    Loads all fasta files from a list or iterable
    """
    return itertools.chain(
        *(
            load_fasta(file_handle)
            for file_handle in files
        )

    )


def load_fasta_rename(file_handle, name_func=None):
    """
    .. versionadded:: 0.3.1

    Renames the header of the sequences using *name_func*, which is called on
    each header. By default, the behaviour is to keep the header to the left of
    the first space (BLAST behaviour).
    """
    if name_func is None:
        def name_func(x): return x.split(' ')[0]

    for seq_id, seq in load_fasta(file_handle):
        yield name_func(seq_id), seq


def load_fasta_prodigal(file_handle):
    """
    .. versionadded:: 0.3.1

    Reads a Prodigal aminoacid fasta file and yields a dictionary with
    basic information about the sequences.

    Arguments:
        file_handle (str, file): passed to :func:`load_fasta`

    Yields:
        dict: dictionary with the information contained in the header, the last
        of the attributes put into key *attr*, while the rest are transformed
        to other keys: seq_id, seq, start, end (genomic), strand, ordinal of
    """

    for seq_id, seq in load_fasta(file_handle):
        prod_seq_id, start, end, strand, attr = seq_id.split(' # ')
        seq_id, idx = prod_seq_id.rsplit('_', 1)

        yield dict(
            prod_seq_id=prod_seq_id,
            seq_id=seq_id,
            seq=seq,
            start=int(start),
            end=int(end),
            strand='+' if strand == '1' else '-',
            idx=int(idx),
            attr=attr
        )


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
    file_handle = mgkit.io.open_file(file_handle, write_mode)

    if wrap is not None:
        seq = '\n'.join(
            seq[pos:pos + wrap]
            for pos in range(0, len(seq), wrap)
        )

    file_handle.write(">{0}\n{1}\n".format(name, seq).encode('ascii'))


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

    def write_func(h, r): return write_fasta_sequence(h, *r)

    mgkit.io.split_write(records, name_mask, write_func, num_files=num_files)
