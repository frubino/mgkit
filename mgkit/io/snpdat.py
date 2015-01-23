"""
SNPDat reader
"""
import logging
import mgkit.utils.sequence

LOG = logging.getLogger(__name__)


class SNPDatRow(object):
    """
    Class containing information ouputted by SNPDat in its result file. One
    instance contains information about a row in the file.

    Attributes:
        chr_name (str): the queried SNPs chromosome ID
        chr_pos (int): queried SNPs genomic location
        in_feat (bool): Whether or not the queried SNP was within a feature
        region (str): Region containing the SNP; either exonic, intronic or
            intergenic
        feat_dist (int): Distance to nearest feature
        feature (str): Either the closest feature to the SNP or the feature
            containing the SNP
        num_features (int): number of different features that the SNP is
            annotated to
        feat_num (int): number of annotations of the current feature
        feat_start (int): Start of feature (bp)
        feat_end (int): End of feature (bp)
        gene_id (str): gene ID for the current feature
        gene_name (str): gene name for the current feature
        transcript_id (str): transcript ID for the current feature
        transcript_name (str): transcript name for the current feature
        exon (tuple): exon that contains the current feature and the total
            number of annotated exons for the gene containing the feature
        strand (str): strand sense of the feature
        ann_frame (str): annotated reading frame (when contained in the GTF)
        frame (str): reading frame estimated by SNPdat
        num_stops (int): estimated number of stop codons in the estimated
            reading frame
        codon (tuple): codon containing the SNP, position in the codon and
            reference base and mutation
        nuc_change (tuple): amino acid for the reference codon and new
            amino acid with the mutation in place
        nuc_ref (str, None): reference nucleotide
        aa_change (str): amino acid for the reference codon and new amino
            acid with the mutation in place
        synonymous (bool): Whether or not the mutation is synonymous
        protein_id (str): protein ID for the current feature
        messages (str): messages in the SNPDat line

    """
    __slots__ = (
        'chr_name',
        'chr_pos',
        'in_feat',
        'region',
        'feat_dist',
        'feature',
        'num_features',
        'feat_num',
        'feat_start',
        'feat_end',
        'gene_id',
        'gene_name',
        'transcript_id',
        'transcript_name',
        'exon',
        'strand',
        'ann_frame',
        'frame',
        'num_stops',
        'codon',
        'nuc_change',
        'nuc_ref',
        'aa_change',
        'synonymous',
        'protein_id',
        'messages'
    )

    def __init__(self, line=None, rev_comp=None):
        if rev_comp is None:
            rev_comp = mgkit.utils.sequence.REV_COMP
        if line:
            line = line.rstrip('\n').split('\t')
            self.chr_name = line[0].lower()  # The queried SNPs chromosome ID
            self.chr_pos = int(line[1])  # The queried SNPs genomic location
            # Whether or not the queried SNP was within a feature
            self.in_feat = True if line[2] == 'Y' else False
            # exonic - in a annotation - check overlapping genes
            # intergenic - between two features
            # NA - no annotation for that contig/snp
            # Region containing the SNP; either exonic, intronic or intergenic
            self.region = line[3].lower()
            # Distance to nearest feature
            if line[4] != 'NA':
                self.feat_dist = int(line[4].replace('*^', ''))
            else:
                self.feat_dist = 'NA'
            # Either the closest feature to the SNP or the feature containing
            # the SNP
            self.feature = line[5].lower()
            # The number of different features that the SNP is annotated to
            self.num_features = 0 if line[6] == 'NA' else int(line[6])
            # The number of annotations of the current feature
            if line[6] == 'NA':
                self.feat_num = 0
            else:
                tuple(int(x) for x in line[7][1:-1].split('/'))
            # Start of feature (bp)
            self.feat_start = 0 if line[6] == 'NA' else int(line[8])
            # End of feature (bp)
            self.feat_end = 0 if line[6] == 'NA' else int(line[9])
            self.gene_id = line[10]  # The gene ID for the current feature
            self.gene_name = line[11]  # The gene name for the current feature
            # The transcript ID for the current feature
            self.transcript_id = line[12]
            # The transcript name for the current feature
            self.transcript_name = line[13]
            # The exon that contains the current feature and the total number
            # of annotated exons for the gene containing the feature
            self.exon = tuple(x for x in line[14][1:-1].split('/'))
            # The strand sense of the feature
            self.strand = line[15].replace('*', '')
            # The annotated reading frame (when contained in the GTF)
            self.ann_frame = line[16]
            # The reading frame estimated by SNPdat
            self.frame = line[17].replace('*', '')
            # The estimated number of stop codons in the estimated reading
            # frame
            if line[6] == 'NA':
                self.num_stops = 0
            else:
                self.num_stops = int(line[18].replace('^', ''))
            # The codon containing the SNP, position in the codon and reference
            # base and mutation
            codon = line[19]
            if codon.startswith('['):
                # first base
                codon = (tuple(codon[1:-3].split('/')), codon[-2], codon[-1])
                if self.strand != '-':
                    self.nuc_change = codon[0][1]
                else:
                    self.nuc_change = rev_comp[codon[0][1]]
                if self.strand != '-':
                    self.nuc_ref = codon[0][0]
                else:
                    self.nuc_ref = rev_comp[codon[0][0]]
            elif codon.endswith(']'):
                # third base
                codon = (codon[0], codon[1], tuple(codon[3:-1].split('/')))
                if self.strand != '-':
                    self.nuc_change = codon[2][1]
                else:
                    self.nuc_change = rev_comp[codon[2][1]]
                if self.strand != '-':
                    self.nuc_ref = codon[2][0]
                else:
                    self.nuc_ref = rev_comp[codon[2][0]]
            elif codon == 'NA':
                codon = None
                self.nuc_change = None
                self.nuc_ref = None
            else:
                # second base
                codon = (codon[0], tuple(codon[2:-2].split('/')), codon[-1])
                if self.strand != '-':
                    self.nuc_change = codon[1][1]
                else:
                    self.nuc_change = rev_comp[codon[1][1]]
                if self.strand != '-':
                    self.nuc_ref = codon[1][0]
                else:
                    self.nuc_ref = rev_comp[codon[1][0]]
            self.codon = codon
            # The amino acid for the reference codon and new amino acid with
            # the mutation in place
            if line[20] == 'NA':
                self.aa_change = (None, None)
            else:
                self.aa_change = tuple(x for x in line[20][1:-1].split('/'))
            # Whether or not the mutation is synonymous
            self.synonymous = True if line[21] == 'Y' else False
            # The protein ID for the current feature
            self.protein_id = line[22]
            # not in a result file if no dbSNP ('-d' options) specified
            # self.rs_id = line[23] #The RS identifier for queries that map to
            # known SNPs
            self.messages = line[23]  # Error messages, warnings etc

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.codon is not None:
            codon = [
                x if isinstance(x, str) else '[{0}/{1}]'.format(*x)
                for x in self.codon
            ]
            codon = "{0}{1}{2}".format(*codon)
        else:
            codon = None
        return "{0}:{1} {2} - [{3}/{4}]".format(self.chr_name, self.chr_pos,
                                                codon, *self.aa_change)


def snpdat_reader(f_handle):
    """
    Simple SNPDat reader.

    f_handle: file handle or string for the SNPDat result file

    :return: generator of SNPDatRow instances
    """

    if isinstance(f_handle, str):
        f_handle = open(f_handle, 'r')
    LOG.info("Reading from file %s", f_handle.name)

    f_handle.readline()  # skips header line

    for line in f_handle:
        try:
            yield SNPDatRow(line)
        except ValueError:
            LOG.critical(line)
            LOG.exception("Error reading line")
            raise ValueError
