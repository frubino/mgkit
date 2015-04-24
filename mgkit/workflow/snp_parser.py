"""
This script parses results of SNPs analysis from any tool for SNP calling [#]_
and integrates them into a format that can be later used for other scripts in
the pipeline.

It integrates coverage and expected number of syn/nonsyn change and taxonomy
from a GFF file, SNP data from a VCF file.

.. note::

    The script accept gzipped VCF files

.. [#] GATK pipeline was the one tested at this time

Changes
*******

.. versionchanged:: 0.1.13
    reworked the internals and the classes used, including options `-m` and `-s`

.. versionchanged:: 0.1.16
    reworkked internals and removed SNPDat, syn/nonsyn evaluation is internal

"""

import HTSeq
import logging
import argparse
import cPickle
from . import utils
from ..io import gff, compressed_handle, fasta
from .. import logger
from ..snps.classes import GeneSNP, SNPType

LOG = logging.getLogger(__name__)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='SNPs analysis, requires a vcf file and SNPDat results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-o',
        '--output-file',
        default='snp_data.pickle',
        type=argparse.FileType('w'),
        help='Ouput file'
    )
    parser.add_argument(
        '-q',
        '--min-qual',
        default=30,
        type=int,
        action='store',
        help='Minimum SNP quality (Phred score)'
    )
    parser.add_argument(
        '-f',
        '--min-freq',
        default=0.01,
        type=float,
        action='store',
        help='Minimum allele frequency'
    )
    parser.add_argument(
        '-r',
        '--min-reads',
        default=4,
        type=int,
        action='store',
        help='Minimum number of reads to accept the SNP'
    )
    parser.add_argument(
        '-g',
        '--gff-file',
        required=True,
        type=argparse.FileType('r'),
        action='store',
        help='GFF file with annotations'
    )
    parser.add_argument(
        '-p',
        '--vcf-file',
        required=True,
        type=argparse.FileType('r'),
        action='store',
        help='Merged VCF file'
    )
    parser.add_argument(
        '-a',
        '--reference',
        required=True,
        type=argparse.FileType('r'),
        help='Fasta file with the GFF Reference'
    )
    parser.add_argument(
        '-m',
        '--samples-id',
        action='append',
        required=True,
        type=str,
        help='the ids of the samples used in the analysis'
    )
    parser.add_argument(
        '-c',
        '--cov-suff',
        action='store',
        default='_cov',
        help="Per sample coverage suffix in the GFF"
    )

    utils.add_basic_options(parser)

    return parser


def init_count_set(annotations):
    LOG.info("Init data structures")

    samples = annotations[0].sample_coverage.keys()

    snp_data = dict(
        (sample, {}) for sample in samples
    )

    for annotation in annotations:

        taxon_id = annotation.taxon_id

        uid = annotation.uid

        sample_coverage = annotation.sample_coverage

        for sample in sample_coverage:
            snp_data[sample][uid] = GeneSNP(
                uid=uid,
                gene_id=annotation.gene_id,
                taxon_id=taxon_id,
                exp_syn=annotation.exp_syn,
                exp_nonsyn=annotation.exp_nonsyn,
                coverage=sample_coverage[sample],
            )

    return snp_data


def check_snp_in_set(samples, snp_data, pos, change, annotations, seq):
    """
    Used by :func:`parse_vcf` to check if a SNP is in a SNPDat result file

    :param iterable samples: list of samples that contain the SNP
    :param dict snp_data: dictionary from :func:`init_count_set` with per
        sample SNPs information
    """

    for annotation in annotations:
        if pos not in annotation:
            continue

        if annotation.is_syn(seq, pos, change):
            snp_type = SNPType.syn
        else:
            snp_type = SNPType.nonsyn

        uid = annotation.uid
        rel_pos = annotation.get_relative_pos(pos)

        for sample in samples:
            snp_data[sample][uid].add_snp(rel_pos, change, snp_type=snp_type)


def parse_vcf(vcf_file, snp_data, min_reads, min_af, min_qual, annotations,
              seqs, line_num=100000):
    """
    Parse VCF file counts synonymous and non-synonymous SNPs

    :param file vcf_file: file handle to a VCF file
    :param dict snp_data: dictionary from :func:`init_count_set` with per
        sample SNPs information
    :param int min_reads: minimum number of reads to accept a SNP
    :param float min_af: minimum allele frequency to accept a SNP
    :param int min_qual: minimum quality (Phred score) to accept a SNP
    :param dict annotations: annotations grouped by their reference sequence
    :param dict seqs: reference sequences
    :param int line_num: the interval in number of lines at which progress
        will be printed
    """
    vcf_handle = HTSeq.VCF_Reader(compressed_handle(vcf_file))

    vcf_handle.parse_meta()
    vcf_handle.make_info_dict()

    #total number of SNPs accepted
    count_tot = 0
    #number of SNPs skipped for low depth
    skip_dp = 0
    #number of SNPs skipped for low allele frequency
    skip_af = 0
    #number of SNPs skipped for low quality
    skip_qual = 0

    for vcf_record in vcf_handle:
        # the SNP is a sequence with no annotations
        if vcf_record.chrom not in annotations:
            continue

        if float(vcf_record.qual) < min_qual:
            #low quality SNP
            skip_qual += 1
            continue

        #unpack info records (needed for vcf_record.info to be a dictionary)
        vcf_record.unpack_info(vcf_handle.infodict)

        #controllare perche' questo controllo e' qui
        if not isinstance(vcf_record.info['DP'], int):
            LOG.warning(vcf_record.info['DP'])

        if vcf_record.info['DP'] < min_reads:
            #not enough reads (depth) for the SNP
            skip_dp += 1
            continue

        allele_freqs = vcf_record.info['AF']

        #if the allele frequency is a single value, make it a list, so
        #the iteration below works anyway
        if isinstance(allele_freqs, float):
            allele_freqs = [allele_freqs]

        # alt is the nucleotidic change
        for allele_freq, change in zip(allele_freqs, vcf_record.alt):
            if allele_freq < min_af:
                #the allele frequency for the SNP is too low, it'll be
                #skipped
                skip_af += 1
                continue

            #the samples that contain the SNP is a string separated by '-'
            samples = [
                sample
                for sample in vcf_record.info['set'].split('-')
            ]
            check_snp_in_set(
                samples,
                snp_data,
                vcf_record.pos.start,
                change,
                annotations[vcf_record.chrom],
                seqs[vcf_record.chrom]
            )
            #increase the total number of snps available
            count_tot += 1

        if vcf_handle.line_no % line_num == 0:
            LOG.info(
                "Line %d, SNPs passed %d; skipped for: qual %d, " +
                "depth %d, freq %d",
                vcf_handle.line_no, count_tot, skip_qual, skip_dp, skip_af
            )


def save_data(output_file, snp_data):
    """
    Pickle data structures to the disk.

    :param str output_file: base name for pickle files
    :param dict snp_data: dictionary from :func:`init_count_set` with per
        sample SNPs information
    """

    LOG.info("Saving sample SNPs to %s", output_file)
    cPickle.dump(snp_data, output_file, -1)


def main():
    "Main function"
    options = set_parser().parse_args()

    #configs log and set log level
    logger.config_log(options.verbose)

    seqs = dict(fasta.load_fasta(options.reference))

    # Loads them as list because it's easier to init the data structure
    annotations = list(gff.parse_gff(options.gff_file))

    LOG.debug("Loaded %d annotations", len(annotations))

    if len(annotations[0].sample_coverage) != len(options.samples_id):
        utils.exit_script(
            "Coverage information was not found for all samples", 2
        )

    snp_data = init_count_set(annotations)

    # Group annotations by their reference sequence
    annotations = gff.group_annotations(
        annotations,
        key_func=lambda x: x.seq_id
    )

    parse_vcf(
        options.vcf_file,
        snp_data,
        options.min_reads,
        options.min_freq,
        options.min_qual,
        annotations,
        seqs
    )

    save_data(options.output_file, snp_data)


if __name__ == '__main__':
    main()
