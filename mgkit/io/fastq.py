"""
Fastq utility functions
"""
import re

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
        (?P<index>[ACTGN]+)"""
"New casava header regex"

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

    print min_qual, max_qual

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
    Write a fastq sequence to file. If the *file_handle* is a string, the file
    will be opened using *write_mode*.

    :param file_handle: file handle or string.
    :param str name: header to write for the sequence
    :param str seq: sequence to write
    :param str qual: quality string
    """
    if isinstance(file_handle, str):
        file_handle = open(file_handle, write_mode)

    file_handle.write(
        "@{name}\n{seq}\n+\n{qual}\n".format(
            name=name,
            seq=seq,
            qual=qual
        )
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
