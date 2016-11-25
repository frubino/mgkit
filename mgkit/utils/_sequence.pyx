import itertools
import numpy
import pandas

def get_kmers(seq, int k):
    cdef int index
    for index in range(0, len(seq) - k + 1):
        yield seq[index:index+k]

def sliding_window(seq, int size, int step):
    cdef int index
    for index in range(0, len(seq) - size + 1, step):
        yield seq[index:index+size]

cdef count_kmer(kmers):
    """
    Simple counter of kmers
    """
    counter = {}
    for kmer in kmers:
        if 'N' in kmer:
            continue
        if kmer not in counter:
            counter[kmer] = 1
        else:
            counter[kmer] += 1
    return counter

def combinations(int k_size):
    "Returns the index for the columns"
    return sorted(''.join(x) for x in itertools.product('ACGT', repeat=k_size))

def sequence_signature(seq, int w_size, int k_size, int step):
    kmer_counts = []
    for subseq in sliding_window(seq, w_size, step):
        kmer_counts.append(count_kmer(get_kmers(subseq, k_size)))
    return kmer_counts

def signatures_matrix(seqs, int w_size, int k_size, int step):
    columns = combinations(k_size)
    n_cols = len(columns)
    index = []
    arrays = []
    cdef int idx
    for seq_id, seq in seqs:
        for idx, counter in enumerate(sequence_signature(seq, w_size, k_size, step)):
            array = numpy.fromiter(
                (counter.get(column, 0) for column in columns),
                numpy.int,
                count=n_cols
            )
            arrays.append(array)
            index.append((seq_id, idx))

    return pandas.DataFrame(numpy.array(arrays), index=pandas.MultiIndex.from_tuples(index), columns=columns, dtype=numpy.int)
