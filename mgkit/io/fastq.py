"""
Fastq utility functions
"""
import re
import numpy
import logging

from .utils import open_file, compressed_handle

LOG = logging.getLogger(__name__)

CASAVA_HEADER_OLD = r"""(?P<machine>\w+-\w+):
        (?P<lane>\d):
        (?P<tile>\d+):
        (?P<xcoord>\d+):
        (?P<ycoord>\d+)\#
        (?P<index>(\d|[ACTGN]{6}))/
        (?P<mate>(1|2))"""
"Old casava header regex"

CASAVA_HEADER_NEW = r"""(?P<machine>[\w-]+):
        (?P<runid>\d+):
        (?P<cellid>\w+):
        (?P<lane>\d):
        (?P<tile>\d+):
        (?P<xcoord>\d+):
        (?P<ycoord>\d+)
        [_ ](?P<mate>\d): # underscore for data from from www.ebi.ac.uk/ena/
        (?P<filter>[YN]):
        (?P<bits>\d+):
        (?P<index>[ACTGN+]+)"""
"New casava header regex, including indices for both forward and reverse"

CASAVA_KHMER = r"""
    (?P<machine>[\w-]+):
    (?P<runid>\d+):
    (?P<cellid>\w+):
    (?P<lane>\d):
    (?P<tile>\d+):
    (?P<xcoord>\d+):
    (?P<ycoord>\d+)
    /(?P<mate>\d)
"""


def check_fastq_type(qualities):
    """
    Trys to guess the type of quality string used in a Fastq file

    :param str qualities: string with the quality scores as in the Fastq file
    :return str: a string with the guessed quality score

    .. note::

        Possible values are the following, classified but the values usually
        used in other softwares:

        * ASCII33: sanger, illumina-1.8
        * ASCII64: illumina-1.3, illumina-1.5, solexa-old

    """
    qualities = [ord(char) for char in qualities]

    max_qual = max(qualities)
    min_qual = min(qualities)

    if (min_qual >= 33) and (max_qual <= 73):
        return 'sanger'
    elif (min_qual >= 33) and (max_qual <= 74):
        return 'illumina-1.8'
    elif (min_qual >= 64) and (max_qual <= 104):
        return 'illumina-1.5'
    elif (min_qual >= 67) and (max_qual <= 104):
        return 'illumina-1.3'
    elif (min_qual >= 59) and (max_qual <= 104):
        return 'solexa-old'


def convert_seqid_to_new(seq_id):
    """
    Convert old seq_id format for Illumina reads to the new found in Casava
    1.8+

    :param str seq_id: seq_id of the sequence (stripped of '@')
    :return str: the new format seq_id

    .. note::

        Example from Wikipedia::

            old casava seq_id:
            @HWUSI-EAS100R:6:73:941:1973#0/1
            new casava seq_id:
            @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCAC

    """
    machine_id, cell_lane, tile_num, x_coord, other = seq_id.split(':')
    y_coord, other = other.split('#')
    idx_multiplex, mate = other.split('/')

    return "{0}:0:FC706VJ:{1}:{2}:{3}:{4} {5}:N:0:{6}".format(
        machine_id,
        cell_lane,
        tile_num,
        x_coord,
        y_coord,
        mate,
        idx_multiplex
    )


def convert_seqid_to_old(seq_id, index_as_seq=True):
    """

    .. deprecated:: 0.3.3

    Convert old seq_id format for Illumina reads to the new found in Casava
    until 1.8, which marks the new format.

    :param str seq_id: seq_id of the sequence (stripped of '@')
    :param bool index_as_seq: if True, the index for the multiplex we'll be
        the sequence found at the end of the new format seq_id. Otherwise, 0
        we'll be used
    :return str: the new format seq_id
    """
    left_part, right_part = seq_id.split(' ')
    machine_id, run_id, cell_id, cell_lane, tile_num, x_coord, \
        y_coord = left_part.split(':')

    mate, fail_filter, control_bits, index_seq = right_part.split(':')

    return "HWUSI-{0}:{1}:{2}:{3}:{4}#{5}/{6}".format(
        machine_id,
        cell_lane,
        tile_num,
        x_coord,
        y_coord,
        index_seq if index_as_seq else '0',
        mate
    )


def write_fastq_sequence(file_handle, name, seq, qual, write_mode='a'):
    """

    .. versionchanged:: 0.3.3
        if *qual* is not a string it's converted to chars (phred33)

    Write a fastq sequence to file. If the *file_handle* is a string, the file
    will be opened using *write_mode*.

    :param file_handle: file handle or string.
    :param str name: header to write for the sequence
    :param str seq: sequence to write
    :param str qual: quality string
    """
    if isinstance(file_handle, str):
        file_handle = open(file_handle, write_mode)

    if isinstance(qual[0], (int, numpy.integer)):
        qual = ''.join(chr(q + 33) for q in qual)

    file_handle.write(
        "@{name}\n{seq}\n+\n{qual}\n".format(
            name=name,
            seq=seq,
            qual=qual
        ).encode('ascii')
    )


def choose_header_type(seq_id):
    """
    Return the guessed compiled regular expression
    :param str seq_id: sequence header to test

    :return: compiled regular expression object or None if no match found
    """
    if re.search(CASAVA_HEADER_OLD, seq_id, re.X) is not None:
        return re.compile(CASAVA_HEADER_OLD, re.X)
    if re.search(CASAVA_HEADER_NEW, seq_id, re.X) is not None:
        return re.compile(CASAVA_HEADER_NEW, re.X)
    if re.search(CASAVA_KHMER, seq_id, re.X) is not None:
        return re.compile(CASAVA_KHMER, re.X)
    return None


def load_fastq(file_handle, num_qual=False):
    """
    .. versionadded:: 0.3.1

    Loads a fastq file and returns a generator of tuples in which the first
    element is the name of the sequence, the second the sequence and the third
    the quality scores (converted in a numpy array if *num_qual* is True).

    .. note::

        this is a simple parser that assumes each sequence is on 4 lines,
        1st and 3rd for the headers, 2nd for the sequence and 4th the quality
        scores

    Arguments:
        file_handle (str, file): fastq file to open, can be a file name or a
            file handle

    Yields:
        tuple: first element is the sequence name/header, the second element is
        the sequence, the third is the quality score. The quality scores are
        kept as a string if *num_qual* is False (default) and converted to a
        numpy array with correct values (0-41) if *num_qual* is True

    Raises:
        ValueError: if the headers in both sequence and quality scores are not
        valid. This implies that the sequence/qualities have carriage returns
        or the file is truncated.

        TypeError: if the qualities are in a format different than sanger
        (min 0, max 40) or illumina-1.8 (0, 41)
    """

    if isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'r')
    else:
        file_handle = compressed_handle(file_handle)

    if getattr(file_handle, 'name', None) is not None:
        LOG.info("Reading fastq file (%s)", file_handle.name)

    check_qual = True

    sequence_count = 0

    while True:
        header1 = file_handle.readline().decode('ascii').strip()
        # Reached the end of the file
        if not header1:
            break

        seq = file_handle.readline().decode('ascii').strip()

        header2 = file_handle.readline().decode('ascii').strip()
        qualities = file_handle.readline().decode('ascii').strip()

        if (header1[0] != '@') or (header2[0] != '+'):
            raise ValueError(
                "The sequence and quality headers are not valid '{}' != '{}'".format(header1, header2)
            )

        header1 = header1[1:]
        header2 = header2[1:]

        if check_qual:
            qual_type = check_fastq_type(qualities)
            if qual_type not in ('sanger', 'illumina-1.8'):
                raise TypeError(
                    'Quality scores are in "{}" format. You must convert them to Sanger (phred33)'.format(qual_type)
                )
            check_qual = False

        if num_qual:
            qualities = numpy.fromiter((ord(x) for x in qualities), numpy.int) - 33

        sequence_count += 1

        yield header1, seq, qualities

    LOG.info("Read %d fastq sequences", sequence_count)


def load_fastq_rename(file_handle, num_qual=False, name_func=None):
    """
    .. versionadded:: 0.3.3

    Mirrors the same functionality in :func:`mgkit.io.fasta.load_fasta_rename`.
    Renames the header of the sequences using *name_func*, which is called on
    each header. By default, the behaviour is to keep the header to the left of
    the first space (BLAST behaviour).
    """

    if name_func is None:
        name_func = lambda header: header.split(' ')[0]

    for header, seq, qual in load_fastq(file_handle, num_qual=num_qual):
        yield name_func(header), seq, qual
