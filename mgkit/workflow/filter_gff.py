"""
Filters GFF annotations in different ways.

Value Filtering
***************

Enables filtering of GFF annotations based on the the values of attributes of a
GFF annotation. The filters are based on equality of numbers (internally
converted into float) and strings, a string contained in the value of an attribute
less or greater than are included as well. The length of annotation has the
attribute *length* and can be tested.

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
`-c` argument to pass a string that defines the function. The string passed must
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

Per Sequence Values
*******************

The *sequence* command allows to filter on a per sequence basis, using
functions such as the median, quantile and mean on attributes like evalue,
bitscore and identity. The file can be passed as sorted already, saving memory
(like in the *overlap* command), but it's not needed to sort the file by strand,
only by the first column.

Coverage Filtering
******************

The *cov* command calculates the coverage of annotations as a measure of the
percentage of each reference sequence length. A minimum coverage percentage can
be used to keep the annotations of sequences that have a greater or equal
coverage than the specified one.

Changes
*******

.. versionadded:: 0.1.12

.. versionchanged:: 0.1.13
    added *--sorted* option

.. versionchanged:: 0.2.0
    changed option *-c* to accept a string to filter overlap

.. versionchanged:: 0.2.5
    added *sequence* command

.. versionchanged:: 0.2.6
    added *length* as attribute and *min*/*max*, and *ge* is the default
    comparison for command *sequence*, *--sort-attr* to *overlap*

.. versionchanged:: 0.3.1
    added *--num-gt* and *--num-lt* to *values* command, added *cov* command

.. versionchanged:: 0.3.4
    moved to use *click* for argument parsing reworked the *values*, *sequence*
    commands

"""
from __future__ import division
from future.utils import viewvalues
import logging
import functools
import click
import pandas
from tqdm import tqdm
import mgkit
from . import utils
from ..io import gff, fasta
from ..filter import gff as filter_gff
from ..utils.common import ranges_length

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@utils.cite_option
def main():
    "Main function"
    pass


def perseq_calc_threshold(annotations, attribute, function, func_arg=None):

    values = pandas.Series(
        annotation.get_attr(attribute, float)
        for annotation in annotations
    )

    if function == 'mean':
        thres = values.mean()
    elif function == 'median':
        thres = values.median()
    elif function == 'quantile':
        thres = values.quantile(func_arg)
    elif function == 'std':
        thres = values.mean() + (func_arg * values.std())
    elif function == 'max':
        thres = values.max()
    elif function == 'min':
        thres = values.min()

    LOG.debug(
        "Threshold for contig %s found: %.2f, using funcion %s and its " +
        "arg %.2f on %d ann.",
        annotations[0],
        thres,
        function,
        func_arg if function in ('quantile', 'std') else 0,
        len(annotations)
    )
    return thres


def find_comparison(comparison):
    if comparison == 'gt':
        return lambda x, y: x > y
    elif comparison == 'ge':
        return lambda x, y: x >= y
    elif comparison == 'lt':
        return lambda x, y: x < y
    elif comparison == 'le':
        return lambda x, y: x <= y


@main.command('sequence', help='Filter on a per sequence basis')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-t', '--sorted', default=False, is_flag=True,
              help='''If the GFF file is sorted (all of a sequence annotations
              are contiguos) can use less memory, `sort -s -k 1,1` can be used''',)
@click.option('-a', '--attribute', default='bitscore', show_default=True,
              type=click.Choice(['evalue', 'bitscore', 'identity', 'length']),
              help='Attribute on which to apply the filter')
@click.option('-f', '--function', show_default=True,
              type=click.Choice(['mean', 'median', 'quantile', 'std', 'max', 'min']),
              default='mean', help='Function for filtering')
@click.option('-l', '--value', default=None, type=click.FLOAT,
              help='Value for the function (used for *std* and *quantile*)')
@click.option('-c', '--comparison', default='ge', show_default=True,
              type=click.Choice(['gt', 'ge', 'lt', 'le']),
              help='Type of comparison (e.g. ge -> greater than or equal to)')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def filter_perseq(verbose, sorted, attribute, function, value, comparison,
                  progress, input_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    if (function in ['std', 'quantile']) and value is None:
        utils.exit_script('When *std* or *quantile* are used, --value must be passed', 1)

    comparison = find_comparison(comparison)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    file_iterator = gff.parse_gff(input_file, gff_type=gff.from_gff)

    if sorted:
        LOG.info("Input GFF is assumed sorted")
        grouped = gff.group_annotations_sorted(
            file_iterator, lambda x: x.seq_id
        )
    else:
        grouped = viewvalues(gff.group_annotations(
            file_iterator, lambda x: x.seq_id
        ))

    if progress:
        grouped = tqdm(grouped)

    for annotations in grouped:
        threshold = perseq_calc_threshold(
            annotations,
            attribute,
            function,
            value
        )

        for annotation in annotations:
            if comparison(annotation.get_attr(attribute, float), threshold):
                annotation.to_file(output_file)


@main.command('cov', help='Filter on a per coverage basis')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-f', '--reference', type=click.File('rb'), required=True,
              help='Reference FASTA file for the GFF')
@click.option('-s', '--strand-specific', default=False, is_flag=True,
              help='If the coverage must be calculated on each strand')
@click.option('-t', '--sorted', default=False, is_flag=True,
              help='Assumes the GFF to be correctly sorted')
@click.option('-c', '--min-coverage', default=0., type=click.FLOAT,
              help='Minimum coverage for the contig/strand')
@click.option('-r', '--rename', default=False, is_flag=True,
              help='Emulates BLAST in reading the FASTA file (keeps only the header before the first space)')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def coverage_command(verbose, reference, strand_specific, sorted, min_coverage,
                     rename, progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)
    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    file_iterator = gff.parse_gff(input_file, gff_type=gff.from_gff)

    if strand_specific:
        def key_func(x): return (x.seq_id, x.strand)
    else:
        def key_func(x): return x.seq_id

    if sorted:
        LOG.info("Input GFF is assumed sorted")
        grouped = gff.group_annotations_sorted(
            file_iterator, key_func
        )
    else:
        grouped = viewvalues(gff.group_annotations(
            file_iterator, key_func
        ))

    seq_iter = fasta.load_fasta_rename(reference) if rename else fasta.load_fasta(reference)

    sequences = {
        seq_id: len(seq)
        for seq_id, seq in seq_iter
    }

    if progress:
        grouped = tqdm(grouped)

    for annotations in grouped:
        covered = ranges_length(
            gff.elongate_annotations(annotations)
        ) / sequences[annotations[0].seq_id] * 100

        if covered < min_coverage:
            continue

        for annotation in annotations:
            annotation.to_file(output_file)


def filter_eq(annotation, attr=None, value=None, conv=None):
    try:
        ann_value = annotation.get_attr(attr, conv=conv)
    except gff.AttributeNotFound:
        return False
    return ann_value == value


def filter_in(annotation, attr=None, value=None, conv=None):
    try:
        ann_value = annotation.get_attr(attr, conv=conv)
    except gff.AttributeNotFound:
        return False
    return value in ann_value


def filter_gt(annotation, attr=None, value=None, conv=None, equal=None):
    try:
        ann_value = annotation.get_attr(attr, conv=conv)
    except gff.AttributeNotFound:
        return False
    if equal:
        return ann_value >= value
    else:
        return ann_value > value


def filter_lt(annotation, attr=None, value=None, conv=None, equal=None):
    try:
        ann_value = annotation.get_attr(attr, conv=conv)
    except gff.AttributeNotFound:
        return False
    if equal:
        return ann_value <= value
    else:
        return ann_value < value


def setup_filters(str_eq, str_in, num_eq, num_ge, num_le, num_gt, num_lt):
    filters = []

    for attr, value in str_eq:
        LOG.info('Filter string (%s = %s)', attr, value)
        filters.append(
            functools.partial(filter_eq, attr=attr, value=value, conv=str)
        )
    for attr, value in str_in:
        LOG.info('Filter string (%s in %s)', value, attr)
        filters.append(
            functools.partial(filter_in, attr=attr, value=value, conv=str)
        )
    for attr, value in num_eq:
        LOG.info('Filter number (%s = %s)', attr, value)
        filters.append(
            functools.partial(filter_eq, attr=attr, value=value, conv=float)
        )
    for attr, value in num_ge:
        LOG.info('Filter number (%s >= %s)', attr, value)
        filters.append(
            functools.partial(filter_gt, attr=attr, value=value, conv=float, equal=True)
        )
    for attr, value in num_gt:
        LOG.info('Filter number (%s > %s)', attr, value)
        filters.append(
            functools.partial(filter_gt, attr=attr, value=value, conv=float, equal=False)
        )
    for attr, value in num_le:
        LOG.info('Filter number (%s <= %s)', attr, value)
        filters.append(
            functools.partial(filter_lt, attr=attr, value=value, conv=float, equal=True)
        )
    for attr, value in num_lt:
        LOG.info('Filter number (%s < %s)', attr, value)
        filters.append(
            functools.partial(filter_lt, attr=attr, value=value, conv=float, equal=False)
        )
    return filters


def validate_params(ctx, param, values, convert=str):
    new_values = []
    for value in values:
        value = value.split(':')
        if len(value) != 2:
            raise click.BadParameter(
                "Wrong key/value format, must be 'key:value' 'key:value'"
            )
        new_values.append((value[0], convert(value[1])))
    return new_values


@main.command('values', help="""Filter based on values""")
@click.option('-v', '--verbose', is_flag=True)
@click.option('--str-eq', multiple=True,
              callback=functools.partial(validate_params, convert=str),
              help='''filter by custom key:value, if the argument is 'key:value' the annotation is kept if it contains an attribute 'key' whose value is exactly 'value' as a string''')
@click.option('--str-in', multiple=True,
              callback=functools.partial(validate_params, convert=str),
              help="Same as '--str-eq' but 'value' is contained in the attribute")
@click.option('--num-eq', multiple=True,
              callback=functools.partial(validate_params, convert=float),
              help="Same as '--str-eq' but 'value' is a number which is equal or greater than")
@click.option('--num-ge', multiple=True,
              callback=functools.partial(validate_params, convert=float),
              help="Same as '--str-eq' but 'value' is a number which is equal or greater than")
@click.option('--num-le', multiple=True,
              callback=functools.partial(validate_params, convert=float),
              help="Same as '--num-ge' but 'value' is a number which is equal or less than")
@click.option('--num-gt', multiple=True,
              callback=functools.partial(validate_params, convert=float),
              help="Same as '--str-eq' but 'value' is a number which is greater than")
@click.option('--num-lt', multiple=True,
              callback=functools.partial(validate_params, convert=float),
              help="Same as '--num-ge' but 'value' is a number which is less than")
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def filter_values(verbose, str_eq, str_in, num_eq, num_ge, num_le, num_gt, num_lt,
                  progress, input_file, output_file):
    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    filters = setup_filters(str_eq, str_in, num_eq, num_ge, num_le, num_gt, num_lt)

    iterator = gff.parse_gff(input_file, gff_type=gff.from_gff)
    if progress:
        iterator = tqdm(iterator)

    for annotation in iterator:
        if all(filter_func(annotation) for filter_func in filters):
            annotation.to_file(output_file)


def make_choose_func(values):
    """Builds the function used to choose between two annotations."""
    values = values.strip()

    LOG.info("Filter function used for overlaps %s", values)

    if values.startswith('+'):
        function = max
        values = values[1:]
    else:
        function = min

    attributes = values.split(',')

    def choose_func(a1, a2):
        return function(
            a1,
            a2,
            key=lambda el: tuple(
                getattr(el, attribute, None) if hasattr(el, attribute) else el.get_attr(attribute, float)
                for attribute in attributes
            )
        )

    return choose_func


@main.command('overlap', help='Use overlapping filter')
@click.option('-v', '--verbose', is_flag=True)
@click.option('-s', '--size', type=click.INT, default=100, show_default=True,
              help='Size of the overlap that triggers the filter')
@click.option('-t', '--sorted', is_flag=True, default=False,
              help='''If the GFF file is sorted (all of a sequence annotations are contiguos and sorted by strand) can use less memory, `sort -s -k 1,1 -k 7,7` can be used''')
@click.option('-c', '--choose-func', default='dbq,bitscore,length',
              help='Function to choose between two overlapping annotations')
@click.option('-a', '--sort-attr', default='bitscore', show_default=True,
              type=click.Choice(['bitscore', 'identity', 'length']),
              help='Attribute to sort annotations before filtering (default bitscore)')
@click.option('--progress', default=False, is_flag=True,
              help="Shows Progress Bar")
@click.argument('input-file', type=click.File('rb'), default='-')
@click.argument('output-file', type=click.File('wb'), default='-')
def filter_overlaps(verbose, size, sorted, choose_func, sort_attr, progress,
                    input_file, output_file):

    mgkit.logger.config_log(level=logging.DEBUG if verbose else logging.INFO)

    choose_func = make_choose_func(choose_func)

    LOG.info(
        'Writing to file (%s)',
        getattr(output_file, 'name', repr(output_file))
    )

    file_iterator = gff.parse_gff(input_file, gff_type=gff.from_gff)

    if sorted:
        LOG.info("Input GFF is assumed sorted")
        grouped = gff.group_annotations_sorted(file_iterator)
    else:
        grouped = viewvalues(gff.group_annotations(file_iterator))

    choose_func = functools.partial(
        filter_gff.choose_annotation,
        overlap=size,
        choose_func=choose_func
    )

    if progress:
        grouped = tqdm(grouped)

    for annotations in grouped:
        filtered = filter_gff.filter_annotations(
            annotations,
            choose_func=choose_func,
            sort_func=lambda x: x.get_attr(sort_attr, float),
            reverse=True
        )
        for annotation in filtered:
            annotation.to_file(output_file)
