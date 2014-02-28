#!/usr/bin/env python
"""
This script parses results of SNPs analysis from any tool for SNP calling [#]_
and `SNPDat <http://code.google.com/p/snpdat/>`_ and integrates them into a
format that can be later used for other scripts in the pipeline.

The data written to disk is two files; each one contains a dictionary with
instances of :class:`snps.GeneSyn` class. One file with suffix '_tot' and one
with suffix '_set':

    * '_tot': overall counts of synonymous/non-synonymous genes from the input
        data
    * '_set': counts specific to each sample

.. note::

    The script accept gzipped VCF files but not SNPDat result files right now

.. [#] GATK pipeline was the one tested at this time
"""

import HTSeq
import gzip
import logging
import argparse
import sys
import cPickle
from . import utils
from ..io import gff
from .. import logger
from ..snps import GeneSyn
from ..io.snpdat import snpdat_reader
from ..taxon import UniprotTaxonomy, MISPELLED_TAXA

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
        '-z',
        '--zipped',
        action='store_true',
        default=False,
        help='If used, the VCF and SNPDat files are expected to be gzipped'
    )
    parser.add_argument(
        '-s',
        '--snpdat-file',
        required=True,
        type=argparse.FileType('r'),
        nargs='+',
        action='store',
        help='SNPDat result file(s)'
    )
    parser.add_argument(
        '-m',
        '--samples-id',
        action='store',
        required=True,
        type=lambda sample: sample.lower(),
        nargs='+',
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
        '-b',
        '--blast',
        action='store_true',
        default=False,
        help="Use BLAST data if available in the GFF file. The taxon used is " +
             "the one found by BLAST."
    )
    parser.add_argument(
        '-t',
        '--taxonomy',
        action='store',
        default=None,
        help="File containing the full Uniprot taxonomy"
    )
    utils.add_basic_options(parser)

    return parser


def reverse_taxon_names(ko2taxon, taxonomy):
    """
    Reverse taxa names

    .. todo::

        take it out and use the taxon_id method of GFF annotations

    """
    LOG.debug("Reversing taxon names for taxa dictionary")
    taxonomy.gen_name_map()
    return dict(
        (ko_idx, taxonomy.find_by_name(taxon_name.replace('#', ' '))[0])
        for ko_idx, taxon_name in ko2taxon.iteritems()
    )


def get_gff_info(gff_file, samples, cov_suff, taxonomy, blast_data):
    """
    Reads the GFF file, reads info about gene expected number of synonymous,
    non-synonymous and per sample coverage.

    Returns:

    * cov_dict: per sample coverage of all genes
    * ko2taxon: dictionary ko_idx->taxon
    * koidx_syn: per gene expected number of synonymous
    * koidx_nonsyn: per gene expected number of non-synonymous

    :param file gff_file: file handle for GFF file
    :param iterable samples: list of sample names
    :param str cov_suff: suffix for the sample coverage attribute in annotations

    :return: cov_dict, ko2taxon, koidx_syn, koidx_nonsyn
    """
    cov_dict = {}

    LOG.info("Getting info from GFF file")

    annotations = gff.load_gff(gff_file)

    ko2taxon = gff.get_attr2attr_map(annotations, keyattr='ko_idx',
                                     valattr='taxon', value_convert=str,
                                     aggr_func=lambda x: x[0])

    ko2taxon_idx = gff.get_attr2attr_map(annotations, keyattr='ko_idx',
                                         valattr='taxon_id', value_convert=int,
                                         aggr_func=lambda x: x[0])

    #correction for the mispelled taxa in old gff file
    #to be taken out after I make sure that all gff files are correct
    ko2taxon_corrected = {}

    for ko_idx, taxon_name in ko2taxon.iteritems():
        try:
            ko2taxon_corrected[ko_idx] = MISPELLED_TAXA[taxon_name]
        except KeyError:
            ko2taxon_corrected[ko_idx] = taxon_name

    ko2taxon = ko2taxon_corrected

    ko2taxon = reverse_taxon_names(ko2taxon, taxonomy)
    ko2taxon.update(ko2taxon_idx)

    if blast_data:
        ko2taxon_blast = gff.get_attr2attr_map(
            annotations,
            keyattr='ko_idx',
            valattr='blast_taxon_idx',
            value_convert=int,
            aggr_func=lambda x: x[0]
        )
        update_ko_taxa(ko2taxon, ko2taxon_blast)

    koidx_syn = gff.get_attr2attr_map(annotations, keyattr='ko_idx',
                                      valattr='exp_syn', value_convert=int,
                                      aggr_func=lambda x: x[0])

    koidx_nonsyn = gff.get_attr2attr_map(annotations, keyattr='ko_idx',
                                         valattr='exp_nonsyn',
                                         value_convert=int,
                                         aggr_func=lambda x: x[0])

    for sample in samples:
        cov_dict[sample] = gff.get_attr2attr_map(annotations, keyattr='ko_idx',
                                                 valattr=sample+cov_suff,
                                                 value_convert=int,
                                                 aggr_func=lambda x: x[0])

    return cov_dict, ko2taxon, koidx_syn, koidx_nonsyn


def init_count_set(samples, ko2taxon, koidx_syn, koidx_nonsyn, cov_dict,
                   taxonomy):
    """
    Init the data structures for the counting of SNPs

    Returns:

    * count_syn: dictionary of :class:`snps.GeneSyn` instances in the form
        ko_idx->GeneSyn
    * count_set: dictionary in the form sample->ko_idx->GeneSyn

    :param iterable samples: list of sample names
    :param dict cov_dict: per sample coverage of all genes
    :param dict ko2taxon: dictionary ko_idx->taxon
    :param dict koidx_syn: per gene expected number of synonymous
    :param dict koidx_nonsyn: per gene expected number of non-synonymous
    :param dict tmap: root taxon map

    :return: count_syn, count_set
    """

    LOG.info("Init data structures")

    count_set = dict(
        (sample, {}) for sample in samples
    )
    count_syn = {}

    #############Add check for KeyError if no coverage was found.
    #or simply use new gff method

    for koidx, taxon_name in ko2taxon.iteritems():

        # print idx, len(ko2taxon)

        taxon_root = taxonomy.get_taxon_root(taxon_name).s_name

        count_syn[koidx] = GeneSyn(
            gid=koidx, taxon=taxon_name,
            exp_syn=koidx_syn[koidx],
            exp_nonsyn=koidx_nonsyn[koidx],
            taxon_root=taxon_root
        )
        for sample in count_set.keys():
            count_set[sample][koidx] = GeneSyn(
                gid=koidx,
                taxon=taxon_name,
                exp_syn=koidx_syn[koidx],
                exp_nonsyn=koidx_nonsyn[koidx],
                coverage=cov_dict[sample][koidx],
                taxon_root=taxon_root
            )

    return count_syn, count_set


def init_count_set2(annotations, taxonomy, prefer_blast):
    """
    Init the data structures for the counting of SNPs

    Returns:

    * count_syn: dictionary of :class:`snps.GeneSyn` instances in the form
        ko_idx->GeneSyn
    * count_set: dictionary in the form sample->ko_idx->GeneSyn

    :param iterable samples: list of sample names
    :param dict cov_dict: per sample coverage of all genes
    :param dict ko2taxon: dictionary ko_idx->taxon
    :param dict koidx_syn: per gene expected number of synonymous
    :param dict koidx_nonsyn: per gene expected number of non-synonymous
    :param dict tmap: root taxon map

    :return: count_syn, count_set
    """

    LOG.info("Init data structures")

    samples = annotations[0].sample_coverage.keys()

    count_set = dict(
        (sample, {}) for sample in samples
    )
    count_syn = {}

    #############Add check for KeyError if no coverage was found.
    #or simply use new gff method

    for annotation in annotations:

        taxon_id = annotation.get_taxon_id(
            taxonomy, prefer_blast=prefer_blast
        )

        taxon_root = taxonomy.get_taxon_root(taxon_id)

        ko_idx = annotation.attributes.ko_idx

        count_syn[ko_idx] = GeneSyn(
            gene_id=ko_idx,
            taxon_id=taxon_id,
            exp_syn=annotation.exp_syn,
            exp_nonsyn=annotation.exp_nonsyn,
            taxon_root=taxon_root
        )

        sample_coverage = annotation.sample_coverage

        for sample in sample_coverage:
            count_set[sample][ko_idx] = GeneSyn(
                gene_id=ko_idx,
                taxon_id=taxon_id,
                exp_syn=annotation.exp_syn,
                exp_nonsyn=annotation.exp_nonsyn,
                coverage=sample_coverage[sample],
                taxon_root=taxon_root
            )

    return count_syn, count_set


def update_ko_taxa(ko2taxon, ko2taxon_blast):
    "Update the taxa dictionary with the blast IDs"
    LOG.info("Updating taxa dictionary with blast id")
    name_dict = {}

    for ko_idx, taxon_id in ko2taxon_blast.iteritems():
        name_dict[ko_idx] = taxon_id

    ko2taxon.update(name_dict)


def load_snpdat_files(snpdat_file, samples, zipped):
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

        if zipped:
            f_handle = gzip.GzipFile(fileobj=f_handle, mode='rb')

        for snp in snpdat_reader(f_handle):
            if snp.nuc_change:
                snp_dat_info[
                    (snp.chr_name, snp.chr_pos, snp.nuc_change, sample)
                ] = snp

    return snp_dat_info


def check_snp_in_set(var_set, snp_dat_info, count_syn, count_set, info_keys):
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
            if snp.synonymous:
                if syn_check is None:
                    count_syn[snp.gene_id].syn += 1
                count_set[var][snp.gene_id].syn += 1
            else:
                if syn_check is None:
                    count_syn[snp.gene_id].nonsyn += 1
                count_set[var][snp.gene_id].nonsyn += 1
            #to make sure we count a snp only once for the total counts
            syn_check = True


def parse_vcf(vcf_file, count_set, count_syn, snp_dat_info, min_reads,
              min_af, min_qual, zipped, line_num=100000):
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
    vcf_handle = HTSeq.VCF_Reader(
        gzip.GzipFile(fileobj=vcf_file, mode='rb') if zipped else vcf_file
    )

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
                (vcf_record.chrom, vcf_record.pos.pos, alt)
            )
            #increase the total number of snps available
            count_tot += 1

        if vcf_handle.line_no % line_num == 0:
            LOG.info("Line %d, SNPs passed %d; skipped for: qual %d, " +
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

    set_fname = "{base}_{type}.{ext}".format(base=base_name, type=set_suff,
                                             ext=ext)
    LOG.info("Saving sample SNPs to %s", set_fname)
    cPickle.dump(count_set, open(set_fname, 'w'), -1)

    tot_fname = "{base}_{type}.{ext}".format(base=base_name, type=tot_suff,
                                             ext=ext)
    LOG.info("Saving total SNPs to %s", tot_fname)
    cPickle.dump(count_syn, open(tot_fname, 'w'), -1)


def main():
    "Main function"
    options = set_parser().parse_args()

    #configs log and set log level
    logger.config_log(options.verbose)

    #the number of SNPDat result files must be the same as the sample ids
    if len(options.snpdat_file) != len(options.samples_id):
        LOG.critical("Number of sample ids and SNPDat result files must be " +
                     "the same")
        sys.exit(1)

    #loads gff file and get all needed information
    #gff_file, samples, cov_suff, taxonomy, blast_data
    taxonomy = UniprotTaxonomy(options.taxonomy)
    # cov_dict, ko2taxon, koidx_syn, koidx_nonsyn = get_gff_info(
    #     options.gff_file,
    #     options.samples_id,
    #     options.cov_suff,
    #     taxonomy,
    #     options.blast
    # )

    annotations = gff.load_gff(options.gff_file)
    LOG.debug(len(annotations))

    count_syn, count_set = init_count_set2(
        annotations,
        taxonomy,
        options.blast
    )

    # del taxonomy  # to free some memory

    snp_dat_info = load_snpdat_files(options.snpdat_file, options.samples_id,
                                     options.zipped)

    parse_vcf(
        options.vcf_file,
        count_set,
        count_syn,
        snp_dat_info,
        options.min_reads,
        options.min_freq,
        options.min_qual,
        options.zipped
    )

    save_data(options.output_file, count_set, count_syn)


if __name__ == '__main__':
    main()
