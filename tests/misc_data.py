import os.path
from nose import SkipTest

from mgkit.io import fasta

data_dir = 'misc_data'
base_dir = os.path.dirname(os.path.abspath(__file__))

### GFF data
gff_name = 'test.gff'

try:
    GFF_FILE = [
        line
        for line in open(os.path.join(base_dir, data_dir, gff_name), 'r')
    ]
except IOError:
    GFF_FILE = None


def setup_gff_data():
    if GFF_FILE is None:
        raise SkipTest(
            'No GFF data found: expecting file "{0}" in dir {1}'.format(
                gff_name,
                data_dir
            )
        )
### HMMER data
hmmer_name = 'test-hmmer-dom.txt'

try:
    HMMER_FILE = [
        line
        for line in open(os.path.join(base_dir, data_dir, hmmer_name), 'r')
    ]
except IOError:
    HMMER_FILE = None


def setup_hmmer_data():
    if HMMER_FILE is None:
        raise SkipTest(
            'No HMMER data found: expecting file "{0}" in dir {1}'.format(
                hmmer_name,
                data_dir
            )
        )

### AA sequences data
aa_name = 'test-seq-aa.fa'

try:
    AA_SEQS = dict(
        fasta.load_fasta(
            os.path.join(base_dir, data_dir, aa_name)
        )
    )
except IOError:
    AA_SEQS = None


def setup_aaseq_data():
    if AA_SEQS is None:
        raise SkipTest(
            'No AA sequences found: expecting file "{0}" in dir {1}'.format(
                aa_name,
                data_dir
            )
        )

### NUC sequences data
nuc_name = 'test-seq-nuc.fa'

try:
    NUC_SEQS = dict(
        fasta.load_fasta(
            os.path.join(base_dir, data_dir, nuc_name)
        )
    )
except IOError:
    NUC_SEQS = None


def setup_nucseq_data():
    if NUC_SEQS is None:
        raise SkipTest(
            'No NUC sequences found: expecting file "{0}" in dir {1}'.format(
                nuc_name,
                data_dir
            )
        )
