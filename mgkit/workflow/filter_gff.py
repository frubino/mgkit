"""
Filters GFF annotations in different ways.

Value Filtering
***************

Enables filtering of GFF annotations based on the the first 8 columns, which
are fixed values as well using the last column which holds information in a
key=value way. There are some predefined key=value filters, like `gene_id`,
but `--str-eq`, `--str-in`, `--num-ge` and `--num-le` allow to make additional
filters.

The functions used to make the filters are located in the module
:mod:`mgkit.filter.gff`, and their names start with `filter_base`, `filter_attr`
and `filter_len`.

.. blockdiag::

    {
        orientation = portrait;

        class mgkit [color = "#e41a1c" , textcolor = 'white', width=200,
        fontsize=15];

        class data [color = "#4daf4a" , textcolor = 'white', width=200,
        fontsize=15];
        "GFF" [class = data, shape = flowchart.input];

        parse_gff [class = "mgkit"];
        setup_filters [class = "mgkit"];
        "Filters" [class = "mgkit", stacked, shape = flowchart.condition];

        "GFF" -> parse_gff -> "Filters";
        setup_filters -> Filters;
        "Filtered Annotations" [class = data, stacked];
        "Filters" -> "Filtered Annotations";

    }

Overlap Filtering
*****************

Filters overlapping annotations using the functions
:func:`mgkit.filter.gff.choose_annotation` and
:func:`mgkit.filter.gff.filter_annotations`, after the annotations are grouped
by both sequence and strand. If the GFF is sorted by sequence name and strand,
the `-t` can be used to make the filtering use less memory. It can be sorted in
Unix using `sort -s -k 1,1 -k 7,7 gff_file`, which applies a stable sort using
the sequence name as the first key and the strand as the second key.

.. note::

    It is also recommended to use::

        export LC_ALL=C

    To speed up the sorting

.. blockdiag::

    {
        orientation = portrait;

        class mgkit [color = "#e41a1c" , textcolor = 'white', width=200,
        fontsize=15];
        class data [color = "#4daf4a" , textcolor = 'white', width=200,
        fontsize=15];
        class software [color = "#377eb8", textcolor = "white", fontsize=15];

        sort [class = software];
        "GFF" [class = data, shape = flowchart.input];
        parse_gff [class = "mgkit"];
        group_annotations [class = "mgkit"];
        filter_annotations [class = "mgkit", shape = flowchart.condition,
        fontsize=12];

        "GFF" -> parse_gff -> sort -> filter_annotations;
        parse_gff -> group_annotations -> filter_annotations;

        "Filtered Annotations" [class = data, stacked];
        filter_annotations -> "Filtered Annotations";
    }

The above digram describes the internals of the script.

The annotations needs first to be grouped by seq_id and strand, forming a group
that can be then be passed to :func:`mgkit.filter.gff.filter_annotations`.
This function:

    #. sort annotations by bit score, from the highest to the lowest
    #. loop over all combination of N=2 annotations:

        #. choose which of the two annotations to discard if they overlap for a
           the required amount of bp (defaults to 100bp)
        #. in which case, the preference is given to the db quality first, than
           the bit score and finally the lenght of annotation, the one with the
           highest values is kept

While the default behaviour is the same, now it is posible to decided the
function used to discard one the two annotations. It is possible to use the
`-c` argument to pass a string that defines the funtion. The string passed must
start with or without a **+**. Using **+** translates into the builtin function
*max* while no **+** translates into *min* from the second character on, any
number of attributes can be used, separated by commas. The attributes, however,
must be one of the properties defined in :class:`mgkit.io.gff.Annotation`,
*bitscore* that returns the value converted in a *float*. Internally the
attributes are stored as strings, so for attributes that have no properties in
the class, such as *evalue*, the `float` builtin is applied.

The tuples built for both annotations are then passed to the comparison
function to be selected and the value returned by it is **discarded**. The
order of the elements in the string is important to define the priority
given to each element in the comparison and the leftmost one has the
highesst priority.

Examples of function strings:

* `-dbq,bitscore,length` becomes max((ann1.dbq, ann1.bitscore, ann1.length),
  (ann2.dbq, ann2.bitscore, ann2.length) - This is default and previously
  only choice
* `-bitscore,length,dbq` uses the same elements but gives lowest priority
  to *dbq*
* `+evalue`: will discard the annotation with the highest *evalue*

Changes
*******

.. versionadded:: 0.1.12

.. versionchanged:: 0.1.13
    added *--sorted* option

.. versionchanged:: 0.2.0
    changed option *-c* to accept a string to filter overlap

"""

import sys
import argparse
import logging
import functools
from .. import logger
from . import utils
from ..io import gff
from ..filter import gff as filter_gff

LOG = logging.getLogger(__name__)


def common_options(parser):
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='Input GFF file, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Output GFF file, defaults to stdout'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_const',
        const=logging.DEBUG,
        default=logging.INFO,
        help='more verbose'
    )


def parse_attr_arg(value, convert=str):
    values = value.split(':')
    if len(values) != 2:
        raise argparse.ArgumentTypeError(
            "Wrong filter format, must be 'key:value' 'key:value'"
        )
    try:
        values[1] = convert(values[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Wrong filter format, cannot convert {0} to number".format(
                values[1]
            )
        )

    return values[0], values[1]


def set_values_parser(parser):
    parser.add_argument(
        '-s',
        '--seq-id',
        action='store',
        type=str,
        help='filter by sequence id'
    )
    parser.add_argument(
        '--strand',
        action='store',
        type=str,
        help='filter by strand',
        choices=['+', '-']
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-sl',
        '--start-lower',
        action='store',
        type=int,
        help='returns only annotations where the start position is less than'
    )
    group.add_argument(
        '-sg',
        '--start-higher',
        action='store',
        type=int,
        help='''returns only annotations where the start position is greater
                than'''
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-el',
        '--end-lower',
        action='store',
        type=int,
        help='''returns only annotations where the end position is equal or
                less than'''
    )
    group.add_argument(
        '-eg',
        '--end-higher',
        action='store',
        type=int,
        help='''returns only annotations where the end position is equal ot
                greater than'''
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-lg',
        '--length',
        action='store',
        type=int,
        help='filter by annotation length equal to or longer than'
    )
    group.add_argument(
        '-ls',
        '--length-short',
        action='store',
        type=int,
        help='filter by annotation length equal to or shorter than'
    )
    parser.add_argument(
        '--source',
        action='store',
        type=str,
        help='filter by source'
    )
    parser.add_argument(
        '-f',
        '--feat-type',
        action='store',
        type=str,
        help='filter by feature type'
    )
    parser.add_argument(
        '-g',
        '--gene-id',
        action='store',
        type=str,
        help='filter by gene_id'
    )
    parser.add_argument(
        '-d',
        '--db',
        action='store',
        type=str,
        help='filter by db'
    )
    parser.add_argument(
        '-q',
        '--db-qual',
        action='store',
        type=int,
        help='filter by db quality equal or greater than'
    )
    parser.add_argument(
        '-b',
        '--bitscore',
        action='store',
        type=float,
        help='filter by bitscore equal or greater than'
    )
    parser.add_argument(
        '-t',
        '--taxon-id',
        action='store',
        type=int,
        help='filter by taxon_id'
    )
    parser.add_argument(
        '--str-eq',
        action='append',
        type=functools.partial(parse_attr_arg, convert=str),
        help='''filter by custom key:value, if the argument is 'key:value' the
             annotation is kept if it contains an attribute 'key' whose value is
             exactly 'value' as a string.
             '''
    )
    parser.add_argument(
        '--str-in',
        action='append',
        type=functools.partial(parse_attr_arg, convert=str),
        help="Same as '--str-eq' but 'value' is contained in the attribute"
    )
    parser.add_argument(
        '--num-ge',
        action='append',
        type=functools.partial(parse_attr_arg, convert=float),
        help="Same as '--str-eq' but 'value' is a number which is equal or greater than"
    )
    parser.add_argument(
        '--num-le',
        action='append',
        type=functools.partial(parse_attr_arg, convert=float),
        help="Same as '--num-ge' but 'value' is a number which is equal or less than"
    )

    common_options(parser)

    parser.set_defaults(func=filter_values)


def set_overlap_parser(parser):
    parser.add_argument(
        '-s',
        '--size',
        action='store',
        type=int,
        help='Size of the overlap that triggers the filter',
        default=100
    )
    parser.add_argument(
        '-t',
        '--sorted',
        action='store_true',
        help='''If the GFF file is sorted (all of a sequence annotations are
                contiguos and sorted by strand) can use less memory,
                `sort -s -k 1,1 -k 7,7` can be used''',
        default=False
    )
    parser.add_argument(
        '-c',
        '--choose-func',
        action='store',
        type=make_choose_func,
        help='Function to choose between two overlapping annotations',
        default='dbq,bitscore,length'
    )

    common_options(parser)

    parser.set_defaults(func=filter_overlaps)


def set_parser():
    """
    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Filter GFF files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()
    parser_b = subparsers.add_parser('values', help='Filter based on values')

    set_values_parser(parser_b)

    parser_o = subparsers.add_parser('overlap', help='Use overlapping filter')

    set_overlap_parser(parser_o)

    utils.add_basic_options(parser)

    return parser


def setup_filters(options):
    filters = []
    #base filters
    base_filters = (
        'taxon_id',
        'seq_id',
        'strand',
        'source',
        'feat_type',
        'gene_id',
        'db'
    )
    for attr in base_filters:
        if getattr(options, attr):
            LOG.info("Filter '%s' = '%s'", attr, getattr(options, attr))
            filters.append(
                functools.partial(
                    filter_gff.filter_base,
                    attr=attr,
                    value=getattr(options, attr)
                )
            )
    #length
    if options.length or options.length_short:

        value = options.length if options.length else options.length_short
        filters.append(
            functools.partial(
                filter_gff.filter_len,
                value=value,
                greater=True if options.length else False
            )
        )
        LOG.info("Filter length %s %s", '>=' if options.length else '<=', value)

    #start position
    if options.start_lower or options.start_higher:
        if options.start_lower:
            value = options.start_lower
        else:
            value = options.start_higher
        filters.append(
            functools.partial(
                filter_gff.filter_base_num,
                attr='start',
                value=value,
                greater=True if options.start_higher else False
            )
        )
        LOG.info(
            "Filter start position %s %s",
            '>=' if options.start_higher else '<=',
            value
        )
    #end position
    if options.end_lower or options.end_higher:
        if options.end_lower:
            value = options.end_lower
        else:
            value = options.end_higher
        filters.append(
            functools.partial(
                filter_gff.filter_base_num,
                attr='end',
                value=value,
                greater=True if options.end_higher else False
            )
        )
        LOG.info(
            "Filter end position %s %s",
            '>=' if options.end_higher else '<=',
            value
        )
    if options.db_qual:
        filters.append(
            functools.partial(
                filter_gff.filter_base_num,
                attr='db_qual',
                value=options.db_qual,
                greater=True
            )
        )
        LOG.info("Filter db quality >= %s", options.db_qual)
    if options.bitscore:
        filters.append(
            functools.partial(
                filter_gff.filter_base_num,
                attr='bitscore',
                value=options.bitscore,
                greater=True
            )
        )
        LOG.info("Filter bitscore >= %s", options.bitscore)
    if options.str_eq:
        for key, value in options.str_eq:
            filters.append(
                functools.partial(
                    filter_gff.filter_attr_str,
                    attr=key,
                    value=value,
                    equal=True
                )
            )
            LOG.info("Filter attribute '%s' = '%s'", attr, value)
    if options.str_in:
        for key, value in options.str_in:
            filters.append(
                functools.partial(
                    filter_gff.filter_attr_str,
                    attr=key,
                    value=value,
                    equal=False
                )
            )
            LOG.info("Filter attribute '%s' contains '%s'", attr, value)
    if options.num_ge:
        for key, value in options.num_ge:
            filters.append(
                functools.partial(
                    filter_gff.filter_attr_num,
                    attr=key,
                    value=value,
                    greater=True
                )
            )
            LOG.info("Filter attribute '%s' >= %s", attr, value)
    if options.num_le:
        for key, value in options.num_le:
            filters.append(
                functools.partial(
                    filter_gff.filter_attr_num,
                    attr=key,
                    value=value,
                    greater=False
                )
            )
            LOG.info("Filter attribute '%s' <= %s", attr, value)

    return filters


def filter_values(options):

    filters = setup_filters(options)

    for annotation in gff.parse_gff(options.input_file, gff_type=gff.from_gff):
        if all(filter(annotation) for filter in filters):
            annotation.to_file(options.output_file)


def make_choose_func(argument):
    """Builds the function used to choose between two annotations."""
    argument = argument.strip()

    LOG.info("Filter function used for overlaps %s", argument)

    if argument.startswith('+'):
        function = max
        argument = argument[1:]
    else:
        function = min

    attributes = argument.split(',')

    choose_func = lambda a1, a2: function(
        a1,
        a2,
        key=lambda el: tuple(
            getattr(el, attribute, None) if hasattr(el, attribute) else el.get_attr(attribute, float)
            for attribute in attributes
        )
    )

    return choose_func


def filter_overlaps(options):

    file_iterator = gff.parse_gff(options.input_file, gff_type=gff.from_gff)

    if options.sorted:
        LOG.info("Input GFF is assumed sorted")
        grouped = gff.group_annotations_sorted(file_iterator)
    else:
        grouped = gff.group_annotations(file_iterator).itervalues()

    choose_func = functools.partial(
        filter_gff.choose_annotation,
        overlap=options.size,
        choose_func=options.choose_func
    )

    for annotations in grouped:
        filtered = filter_gff.filter_annotations(
            annotations,
            choose_func=choose_func,
            sort_func=lambda x: x.bitscore,
            reverse=True
        )
        for annotation in filtered:
            annotation.to_file(options.output_file)


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)

    options.func(options)
