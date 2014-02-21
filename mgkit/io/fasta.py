# coding=utf8
"""
Simple fasta parser and a few utility functions
"""
import collections
# import gff
import logging

TEXT_WRAP_DEFAULT = 80
"""
Number of characters to use for wrapping
"""

LOG = logging.getLogger(__name__)


def load_fasta(f_handle):
    """
    Loads a fasta file and returns a generator of tuples in which the first
    element is the name of the sequence and the second the sequence
    """
    if isinstance(f_handle, str):
        f_handle = open(f_handle, 'r')

    if getattr(f_handle, 'name', None) is not None:
        LOG.info("Reading fasta file %s", f_handle.name)

    cur_name = ""
    last_name = ""
    nseq = 0
    # Better use ''.join(list) than seq+=seq, it's faster
    cur_seq = []
    # main loop to read file's sequences into an Alignment Object
    for line in f_handle:
        if line.startswith('>'):
            if cur_seq != []:
                # start of next sequence
                yield cur_name, ''.join(cur_seq)
                nseq += 1
                cur_seq = []
            # save previous name for loop's else clause
            last_name = cur_name
            #fasta classico, prende tutti i caratteri tranne il primo ">"
            #e il whitespace a destra
            cur_name = line[1:].rstrip()
        else:
            cur_seq.append(line.rstrip())
    else:
        # handles cases where there's no newline after last
        # portion of the sequence at EOF
        if nseq != 0:
            if cur_name != last_name:
                yield cur_name, ''.join(cur_seq)
        # case in which only one sequence is present
        else:
            yield cur_name, ''.join(cur_seq)
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

    file_handle.write(">{0}\n".format(name))

    if wrap is None:
        file_handle.write("{0}\n".format(seq))
        return

    for pos in xrange(0, len(seq), wrap):
        file_handle.write("{0}\n".format(seq[pos:pos+wrap]))


def grep_fasta(f_in, f_out, grep_file):
    f = open(f_in, 'r')
    Sequence = collections.namedtuple('Sequence', "name seq")
    LOG.info("Loading grep file %s", grep_file)
    contigs = set([x.strip() for x in open(grep_file, 'r')])
    LOG.info("Number of contigs in grep file %d", len(contigs))
    seqs = [Sequence(name, seq) for name, seq in load_fasta(f)]
    f = open(f_out, 'w')
    count = 0
    for name, seq in seqs:
        if name in contigs:
            count += 1
            f.write(">{0}\n{1}\n".format(name, seq))
    LOG.info("Number of grepped sequences %d", count)


# def filter_fasta(f_in, f_out, gff_file):
#     f = open(f_in, 'r')
#     Sequence = collections.namedtuple('Sequence', "name seq")

#     LOG.info("Loading GFF annotations %s", gff_file)
#     contigs = set(x.seq_id for x in gff.load_gff(open(gff_file, 'r')))
#     seqs = [Sequence(name, seq) for name, seq in load_fasta(f) if name in contigs]
#     f = open(f_out, 'w')

#     print "Writing %d sequences to file %s" % (len(seqs), f_out)
#     for x in seqs:
#         f.write(">{0}\n{1}\n".format(x.name, x.seq))


# def contigs_stats(fname, all_length, gff_file=None):
#     f = open(fname, 'r')
#     Sequence = collections.namedtuple('Sequence', "name seq")
#     if gff_file is None:
#         seqs = [Sequence(name, seq) for name, seq in load_fasta(f)]
#     else:
#         print "Loading GFF annotations", gff_file
#         contigs = set(x.seq_id for x in gff.load_gff(open(gff_file, 'r')))
#         seqs = [Sequence(name, seq) for name, seq in load_fasta(f) if name in contigs]

#     all_len = [len(x.seq) for x in seqs]
#     avg_len = sum(all_len) / float(len(seqs))
#     print "Number of sequences %d with an average length of %f, minimum %f, maximum %d" % (len(seqs), avg_len, min(all_len), max(all_len))
#     f = open(all_length, 'w')
#     for x in all_len:
#         f.write("{0}\n".format(x))
