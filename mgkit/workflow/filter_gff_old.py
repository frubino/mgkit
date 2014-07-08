"""
Filter a GFF file
"""
import sys
import argparse
import functools
import logging
import pickle
from mgkit.io import gff
from mgkit import logger
from mgkit.filter.gff_old import *
from . import utils


def set_parser():
    "set options"
    parser = argparse.ArgumentParser(
        description='Filter GFF files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    utils.add_basic_options(parser)

    parser.add_argument(
        '-q', '--seq-id', nargs='+', metavar='SEQID',
        help='more than one SEQID can be specified'
    )
    parser.add_argument('-d', '--strand', choices=['+', '-'])

    group = parser.add_argument_group('HMMR options')
    group.add_argument(
        '-s', '--score', metavar='EVALUE', type=float,
        help='Annotation is not saved if its score > EVALUE'
    )
    group.add_argument(
        '-b', '--bit-score', type=float,
        help='Annotation is not saved if its bit-score < BIT_SCORE'
    )
    group.add_argument(
        '-f', '--profile-len', type=argparse.FileType('r'),
        action='append',
        help="""pickle file(s) containing the profile length"""
    )
    group.add_argument(
        '-p', '--gene-percentage', type=float, default=0.4,
        help="""Hit must be of length equal or greater than a percentange of
        the profile length"""
    )

    group = parser.add_argument_group('Kegg options')
    group.add_argument(
        '-t', '--taxon', nargs='+',
        help='more than one TAXON can be specified'
    )
    group.add_argument(
        '-k', '--ko-id', nargs='+', metavar='KO',
        help='more than one KO can be specified'
    )
    group.add_argument(
        '-K', '--ko-idx', nargs='*', metavar='KO_IDX',
        help='more than one KO_IDX can be specified'
    )
    group.add_argument(
        '-r', '--only-reviewed', action='store_true',
        help='include only reviewed profiles'
    )
    # parser.add_argument('-c', '--case-insesitive', action='store_true',
                          # help='all string lookups are case insensitive')
    group.add_argument(
        '-e', '--description',
        help='looks for DESCRIPTION in the annotation (case insensitive)'
    )

    group = parser.add_argument_group('overlapping filter options')
    group.add_argument(
        '-l', '--filter-overlap', action='store_true',
        help='filter overlapping annotations, maintain the one with the best ' +
        'score'
    )
    group.add_argument(
        '-z', '--threshold', type=float, default=0.8,
        help='threshold used to choose between overlapping regions ' +
        '(percentange of average lenght)'
    )
    group.add_argument(
        '-g', '--same-gene', default=False, action='store_true',
        help='specify that the overlaps apply only to the same gene'
    )
    group.add_argument(
        '-n', '--both-strands', default=True, action='store_false',
        help='specify that the overlaps do not have to check the strand'
    )

    group = parser.add_argument_group('file options')
    group.add_argument('-i', '--input-file', nargs='?',
                       type=argparse.FileType('r'),
                       default='-', required=True)
    group.add_argument('-F', '--pattern-file', type=argparse.FileType('r'))
    group.add_argument('-o', '--output-file', nargs='?',
                       type=argparse.FileType('w'),
                       default=sys.stdout)

    return parser


def set_filter_list(options):
    "Prepare filter list"
    filter_list = []

    log = logging.getLogger(__name__)

    if options.score:
        log.info(
            "Only include annotations with score greater than or equal to %f",
            options.score
        )
        filter_list.append(
            functools.partial(filter_by_score, options.score)
        )
    if options.bit_score:
        log.info(
            "Only include annotations with bit-score greater than or " +
            "equal to %f",
            options.bit_score
        )
        filter_list.append(
            functools.partial(filter_by_bit_score, options.bit_score)
        )
    if options.taxon:
        log.info(
            "Only include annotation for taxa: %s",
            ','.join([x.lower() for x in options.taxon])
        )
        filter_list.append(
            functools.partial(
                filter_by_taxon, [x.lower() for x in options.taxon]
            )
        )
    if options.seq_id:
        log.info(
            "Only include annotation for seq_id: %s",
            ','.join([x.lower() for x in options.seq_id])
        )
        filter_list.append(
            functools.partial(
                filter_by_seq_id, [x.lower() for x in options.seq_id]
            )
        )
    if not options.ko_idx is None:
        patterns = options.ko_idx
        if options.pattern_file:
            log.debug("Using pattern file %s", options.pattern_file.name)
            patterns = [x.rstrip() for x in options.pattern_file]
        log.info(
            "Only include annotation for ko_idx: %s...",
            ','.join([x.lower() for x in patterns[:10]])
        )
        filter_list.append(
            functools.partial(filter_by_ko_idx, [x.lower() for x in patterns])
        )
    if options.ko_id:
        log.info(
            "Only include annotation for kos: %s",
            ','.join([x.lower() for x in options.ko_id])
        )
        filter_list.append(
            functools.partial(
                filter_by_ko_id, [x.lower() for x in options.ko_id]
            )
        )
    if options.strand:
        log.info("Only include strand '%s'", options.strand)
        filter_list.append(
            functools.partial(filter_by_strand, options.strand)
        )
    if options.description:
        log.info(
            "Only include annotation whose description contains '%s'",
            options.description
        )
        filter_list.append(
            functools.partial(
                filter_by_description, options.description.lower()
            )
        )
    if options.only_reviewed:
        log.info("Only include reviewed profiles")
        filter_list.append(filter_by_reviewed)
    if options.profile_len:
        log.info(
            "Set filtering with profile length (%.2f) file(s) %s",
            options.gene_percentage,
            ','.join(x.name for x in options.profile_len)
        )
        ko_len = {}
        for dictionary in options.profile_len:
            ko_len.update(pickle.load(dictionary))
        filter_list.append(
            functools.partial(
                filter_by_hit_length,
                ko_len,
                options.gene_percentage
            )
        )

    return filter_list


def main():
    "Main function"
    options = set_parser().parse_args()
    logger.config_log()
    log = logging.getLogger(__name__)
    # print set_parser().parse_args('-h -q con1 con2 -k arc bac -i -'.split())
    filter_list = set_filter_list(options)

    if options.filter_overlap:
        contigs = {}

    for idx, line in enumerate(options.input_file):
        if (idx + 1) % 100000 == 0:
            log.info("Filtered %d lines", idx + 1)
        annotation = gff.GFFKegg(line)
        if all(x(annotation) for x in filter_list):
            if options.filter_overlap:
                try:
                    contigs[annotation.seq_id].append(annotation)
                except KeyError:
                    contigs[annotation.seq_id] = [annotation]
            else:
                options.output_file.write(str(annotation))
    if options.filter_overlap:
        length = len(contigs)
        log.info(
            'Start overlap filtering with options: overlap threshold (%f), ' +
            'check only same gene (%s) check strand (%s)',
            options.threshold,
            options.same_gene,
            options.both_strands
        )
        for idx, annotations in enumerate(contigs.itervalues()):
            log.debug(
                "(%d/%d) % 14s - Num. of features: %d", idx+1, length,
                annotations[0].seq_id, len(annotations)
            )
            filtered = filter_overlapping(
                annotations,
                options.threshold,
                options.same_gene,
                options.both_strands
            )

            for annotation in filtered:
                options.output_file.write(str(annotation))

if __name__ == '__main__':
    main()
