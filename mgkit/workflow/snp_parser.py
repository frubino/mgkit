"""
.. versionchanged:: 0.1.13
    reworked the internals and the classes used

This script parses results of SNPs analysis from any tool for SNP calling [#]_
and `SNPDat <http://code.google.com/p/snpdat/>`_ and integrates them into a
format that can be later used for other scripts in the pipeline.

It integrates coverage and expected number of syn/nonsyn change and taxonomy
from a GFF file, SNP data from a VCF file, syn/nonsyn information from SNPDat.

The data written to disk is two files; each one contains a dictionary with
instances of :class:`snps.GeneSyn` class. One file with suffix '_tot' and one
with suffix '_set':

    * '_tot': overall counts of synonymous/non-synonymous genes from the input
        data
    * '_set': counts specific to each sample

.. note::

    The script accept gzipped VCF and SNPDat result files

.. [#] GATK pipeline was the one tested at this time
"""

import HTSeq
import logging
import argparse
import cPickle
from . import utils
from ..io import gff, compressed_handle
from .. import logger
from ..snps.classes import GeneSNP, SNPType
from ..io.snpdat import snpdat_reader
from ..taxon import UniprotTaxonomy

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
        default='snps_out',
        action='store',
        help='Ouput basename: two file are created, with suffix _tot.pickle ' +
             'and _set.pickle'
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
        '-s',
        '--snpdat-file',
        required=True,
        type=argparse.FileType('r'),
        action='append',
        help='SNPDat result file(s)'
    )
    parser.add_argument(
        '-m',
        '--samples-id',
        action='append',
        required=True,
        type=lambda sample: sample.lower(),
        help='the ids of the samples used in the analysis'
    )
    parser.add_argument(
        '-c',
        '--cov-suff',
        action='store',
        default='_cov',
        help="Per sample coverage suffix in the GFF. (e.g. if '_cov' the " +
             "attribute in the annotation must be 'sample_cov'"
    )
    parser.add_argument(
        '-t',
        '--taxonomy',
        action='store',
        default=None,
        help="File containing the full Uniprot taxonomy"
    )
    parser.add_argument(
        '-d',
        '--corr-old',
        action='store_true',
        default=False,
        help="Correct old annotations taxonomy"
    )
    utils.add_basic_options(parser)

    return parser


def init_count_set(annotations, taxonomy, correct_ann=False):
    LOG.info("Init data structures")

    samples = annotations[0].sample_coverage.keys()

    count_set = dict(
        (sample, {}) for sample in samples
    )
    count_syn = {}

    if correct_ann:
        gff.correct_old_annotations(annotations, taxonomy)

    for annotation in annotations:

        taxon_id = annotation.taxon_id

        ko_idx = annotation.uid

        count_syn[ko_idx] = GeneSNP(
            uid=ko_idx,
            gene_id=annotation.gene_id,
            taxon_id=taxon_id,
            exp_syn=annotation.exp_syn,
            exp_nonsyn=annotation.exp_nonsyn
        )

        sample_coverage = annotation.sample_coverage

        for sample in sample_coverage:
            count_set[sample][ko_idx] = GeneSNP(
                uid=ko_idx,
                gene_id=annotation.gene_id,
                taxon_id=taxon_id,
                exp_syn=annotation.exp_syn,
                exp_nonsyn=annotation.exp_nonsyn,
                coverage=sample_coverage[sample],
            )

    return count_syn, count_set


def load_snpdat_files(snpdat_file, samples):
    """
    Loads information from SNPDat result files

    The returned dictionary uses as keys same attributes from the
    :class:`snps.SNPDatRow` instance, adding the sample name to make it unique
    among samples.

    :param iterable snpdat_file: list of file handles of SNPDat result files
    :param iterble samples: sample names; must be in the same order of the
        snpdat_file

    :return dict: dictionary in the form
        (chr_name, chr_pos, nuc_change, sample)->:class:`snps.SNPDatRow`

    """
    snp_dat_info = {}

    for f_handle, sample in zip(snpdat_file, samples):

        f_handle = compressed_handle(f_handle)

        for snp in snpdat_reader(f_handle):
            if snp.nuc_change:
                snp_dat_info[
                    (snp.chr_name, snp.chr_pos, snp.nuc_change, sample)
                ] = snp

    return snp_dat_info


def check_snp_in_set(var_set, snp_dat_info, count_syn, count_set, info_keys,
                     annotations, start):
    """
    Used by :func:`parse_vcf` to check if a SNP is in a SNPDat result file

    :param iterable var_set: list of samples that contain the SNP
    :param dict snp_dat_info: dictionary with SNPDat result file information
        from :func:`load_snpdat_files`
    :param dict count_set: dictionary from :func:`init_count_set` with per
        sample SNPs information
    :param dict count_syn: dictionary from :func:`init_count_set` with total
        SNPs information
    :param tuple info_keys: tuple with the keys used to retrieve the data from
        'snp_dat_info' to which is added the sample name

    """
    #used to make sure that in the overall count ('count_syn') a SNP is counted
    #once in the al loop
    syn_check = None
    for var in var_set:
        try:
            # Data available in the snpdat results
            snp = snp_dat_info[info_keys + (var,)]
        except KeyError:
            continue
        if snp.region == 'exonic':
            gene_start = start - annotations[snp.gene_id].start + 1
            if snp.synonymous:
                snp_tuple = (gene_start, info_keys[2], SNPType.syn)
                if syn_check is None:
                    count_syn[snp.gene_id].add_snp(*snp_tuple)
                count_set[var][snp.gene_id].add_snp(*snp_tuple)
            else:
                snp_tuple = (gene_start, info_keys[2], SNPType.nonsyn)
                if syn_check is None:
                    count_syn[snp.gene_id].add_snp(*snp_tuple)
                count_set[var][snp.gene_id].add_snp(*snp_tuple)
            #to make sure we count a snp only once for the total counts
            syn_check = True


def parse_vcf(vcf_file, count_set, count_syn, snp_dat_info, min_reads,
              min_af, min_qual, annotations, line_num=100000):
    """
    Parse VCF file counts synonymous and non-synonymous SNPs

    :param file vcf_file: file handle to a VCF file
    :param dict count_set: dictionary from :func:`init_count_set` with per
        sample SNPs information
    :param dict count_syn: dictionary from :func:`init_count_set` with total
        SNPs information
    :param dict snp_dat_info: returned from :func:`load_snpdat_files`
    :param int min_reads: minimum number of reads to accept a SNP
    :param float min_af: minimum allele frequency to accept a SNP
    :param int min_qual: minimum quality (Phred score) to accept a SNP
    :param bool zipped: if True the file is gzipped
    :param int line_num: the interval in number of lines at which progress will
        be printed
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

        for allele_freq, alt in zip(allele_freqs, vcf_record.alt):
            if allele_freq < min_af:
                #the allele frequency for the SNP is too low, it'll be
                #skipped
                skip_af += 1
                continue

            #the samples that contain the SNP is a string separated by '-'
            var_set = [
                sample.lower()
                for sample in vcf_record.info['set'].split('-')
            ]
            check_snp_in_set(
                var_set,
                snp_dat_info,
                count_syn,
                count_set,
                #keep pos.pos for backward compatibility
                (vcf_record.chrom, vcf_record.pos.pos, alt),
                annotations,
                vcf_record.pos.start
            )
            #increase the total number of snps available
            count_tot += 1

        if vcf_handle.line_no % line_num == 0:
            LOG.info(
                "Line %d, SNPs passed %d; skipped for: qual %d, " +
                "depth %d, freq %d; tot syn/nonsyn: %r",
                vcf_handle.line_no, count_tot, skip_qual, skip_dp, skip_af,
                (
                    sum(x.syn for x in count_syn.values()),
                    sum(x.nonsyn for x in count_syn.values())
                )
            )


def save_data(base_name, count_set, count_syn, tot_suff='tot', set_suff='set',
              ext='pickle'):
    """
    Pickle data structures to the disk.

    :param str base_name: base name for pickle files
    :param dict count_set: dictionary from :func:`init_count_set` with per
        sample SNPs information
    :param dict count_syn: dictionary from :func:`init_count_set` with total
        SNPs information
    :param str tot_suff: suffix added to 'base_name' for 'count_syn' file name
    :param str set_suff: suffix added to 'base_name' for 'count_set' file name
    :param str ext: file extension

    """

    set_fname = "{base}_{type}.{ext}".format(
        base=base_name,
        type=set_suff,
        ext=ext
    )
    LOG.info("Saving sample SNPs to %s", set_fname)
    cPickle.dump(count_set, open(set_fname, 'w'), -1)

    tot_fname = "{base}_{type}.{ext}".format(
        base=base_name,
        type=tot_suff,
        ext=ext
    )
    LOG.info("Saving total SNPs to %s", tot_fname)
    cPickle.dump(count_syn, open(tot_fname, 'w'), -1)


def main():
    "Main function"
    options = set_parser().parse_args()

    #configs log and set log level
    logger.config_log(options.verbose)

    #the number of SNPDat result files must be the same as the sample ids
    if len(options.snpdat_file) != len(options.samples_id):
        utils.exit_script(
            "Number of sample ids and SNPDat result files must be the same",
            1
        )

    if options.taxonomy:
        taxonomy = UniprotTaxonomy(options.taxonomy)
    else:
        taxonomy = options.taxonomy
        if options.corr_old:
            utils.exit_script(
                "To correct old annotation the taxonomy is required",
                3
            )

    annotations = list(gff.parse_gff(options.gff_file))
    LOG.debug("Loaded %d annotations", len(annotations))

    if len(annotations[0].sample_coverage) != len(options.samples_id):
        utils.exit_script("Coverage information was not found for all samples", 2)

    count_syn, count_set = init_count_set(
        annotations,
        taxonomy,
        options.corr_old
    )

    snp_dat_info = load_snpdat_files(
        options.snpdat_file, options.samples_id
    )

    annotations = dict((x.uid, x) for x in annotations)

    parse_vcf(
        options.vcf_file,
        count_set,
        count_syn,
        snp_dat_info,
        options.min_reads,
        options.min_freq,
        options.min_qual,
        annotations
    )

    save_data(options.output_file, count_set, count_syn)


if __name__ == '__main__':
    main()
