"""
.. versionadded:: 0.1.13

Misc functions for count data
"""
import sys
import logging
import itertools
import functools
from collections import Counter
from future.utils import viewitems

from mgkit.filter import taxon as tx_filters
from mgkit.io import open_file
import mgkit.simple_cache

import pandas

LOG = logging.getLogger(__name__)

SKIP = set(
    [
        'no_feature',
        'ambiguous',
        'too_low_aQual',
        'not_aligned',
        'alignment_not_unique'
    ]
)


def load_htseq_counts(file_handle, conv_func=int):
    """
    .. versionchanged:: 0.1.15
        added *conv_func* parameter

    Loads an HTSeq-count result file

    Arguments:

        file_handle (file or str): file handle or string with file name
        conv_func (func): function to convert the number from string, defaults
            to *int*, but *float* can be used as well

    Yields:
        tuple: first element is the gene_id and the second is the count

    """

    if (sys.version_info[0] == 2) and isinstance(file_handle, unicode):
        file_handle = open_file(file_handle, 'rb')
    elif isinstance(file_handle, str):
        file_handle = open_file(file_handle, 'rb')

    if getattr(file_handle, 'name', None) is not None:
        LOG.info("Loading HTSeq-count file %s", file_handle.name)

    for line in file_handle:
        line = line.decode('ascii')
        gene_id, count = line.rstrip().split('\t')

        if line.startswith('__') or (gene_id in SKIP):
            continue

        yield gene_id, conv_func(count)


def batch_load_htseq_counts(count_files, samples=None, cut_name=None):
    """
    Loads a list of htseq count result files and returns a DataFrame
    (IDxSAMPLE)

    The sample names are names are the file names if *samples* and *cut_name*
    are *None*, supplying a list of sample names with *samples* is the
    preferred way, and *cut_name* is used for backward compatibility and as an
    option in cases a string replace is enough.

    Arguments:
        count_files (file or str): file handle or string with file name
        samples (iterable): list of sample names, in the same order as
            *count_files*
        cut_name (str): string to delete from the the file names to get the
            sample names

    Returns:
        pandas.DataFrame: with sample names as columns and gene_ids as index
    """
    counts = {}

    iterator = itertools.zip_longest(
        count_files,
        [] if samples is None else samples
    )

    for fname, sample in iterator:
        if sample is None:
            if cut_name is not None:
                sample = fname.replace(cut_name, '')
            else:
                sample = fname

        counts[sample] = pandas.Series(load_htseq_counts(fname))

    return pandas.DataFrame.from_dict(counts)


def get_uid_info(info_dict, uid):
    """
    Simple function to get a value from a dictionary of tuples
    (gene_id, taxon_id)
    """
    return info_dict[uid]


def get_uid_info_ann(annotations, uid):
    """
    Simple function to get a value from a dictionary of annotations
    """
    annotation = annotations[uid]
    return annotation.gene_id, annotation.taxon_id


def filter_counts(counts_iter, info_func, gfilters=None, tfilters=None):
    """
    Returns counts that pass filters for each *uid* associated *gene_id* and
    *taxon_id*.

    Arguments:
        counts_iter (iterable): iterator that yields a tuple (uid, count)
        info_func (func): function accepting a *uid* that returns a tuple
            *(gene_id, taxon_id)*
        gfilters (iterable): list of filters to apply to each *uid* associated
            *gene_id*
        tfilters (iterable): list of filters to apply to each *uid* associated
            *taxon_id*

    Yields:
        tuple: *(uid, count)* that pass filters
    """
    for uid, count in counts_iter:
        try:
            gene_id, taxon_id = info_func(uid)
        except KeyError:
            continue

        if tfilters is not None:
            if taxon_id is None:
                continue

            if not all(tfilter(taxon_id) for tfilter in tfilters):
                continue

        if gfilters is not None:
            if not all(gfilter(taxon_id) for gfilter in gfilters):
                continue

        yield uid, count


def map_counts(counts_iter, info_func, gmapper=None, tmapper=None, index=None,
               uid_used=None):
    """
    .. versionchanged:: 0.1.14
        added *index* parameter

    .. versionchanged:: 0.1.15
        added *uid_used* parameter

    Maps counts according to the gmapper and tmapper functions. Each mapped
    gene ID count is the sum of all uid that have the same ID(s). The same is
    true for the taxa.

    Arguments:
        counts_iter (iterable): iterator that yields a tuple (uid, count)
        info_func (func): function accepting a *uid* that returns a tuple
            *(gene_id, taxon_id)*
        gmapper (func): fucntion that accepts a *gene_id* and returns a list
            of mapped IDs
        tmapper (func): fucntion that accepts a *taxon_id* and returns a new
            *taxon_id*
        index (None, str): if None, the index of the Series if
            *(gene_id, taxon_id)*, if a str, it can be either *gene* or
            *taxon*, to specify a single value
        uid_used (None, dict): an empty dictionary in which to store the *uid*
            that were assigned to each key of the returned pandas.Series. If
            *None*, no information is saved

    Returns:
        pandas.Series: array with MultiIndex *(gene_id, taxon_id)* with the
        mapped counts
    """
    mapped_counts = {}

    for uid, count in counts_iter:
        try:
            gene_id, taxon_id = info_func(uid)
        except KeyError:
            continue

        if gmapper is not None:
            gene_ids = gmapper(gene_id)
        else:
            gene_ids = [gene_id]

        if tmapper is not None:
            taxon_id = tmapper(taxon_id)
            if taxon_id is None:
                continue

        for map_id in gene_ids:
            if index is None:
                key = (map_id, taxon_id)
            elif index == 'gene':
                key = map_id
            elif index == 'taxon':
                key = taxon_id

            try:
                mapped_counts[key] += count
                if uid_used is not None:
                    uid_used[key].add(uid)
            except KeyError:
                mapped_counts[key] = count
                if uid_used is not None:
                    uid_used[key] = set([uid])

    return pandas.Series(mapped_counts)


def map_taxon_id_to_rank(taxonomy, rank, taxon_id, include_higher=True):
    """
    Maps a *taxon_id* to the request taxon rank. Returns *None* if
    *include_higher* is False and the found rank is not the one requested.

    Internally uses :meth:`mgkit.taxon.Taxonomy.get_ranked_taxon`

    Arguments:
        taxonomy: taxonomy instance
        rank (str): taxonomic rank requested
        taxon_id (int): taxon_id to map
        include_higher (bool): if False, any rank different than the requested
            one is discarded

    Returns:
        (int, None): if the mapping is successful, the ranked taxon_id is
        returned, otherwise *None* is returned
    """
    if taxon_id is None:
        return None

    ranked = taxonomy.get_ranked_taxon(taxon_id, rank)

    if (include_higher is False) and (ranked.rank != rank):
        return None
    return ranked.taxon_id


def map_gene_id_to_map(gene_map, gene_id):
    """
    Function that extract a list of gene mappings from a dictionary and returns
    an empty list if the *gene_id* is not found.
    """
    return gene_map.get(gene_id, [])


def load_sample_counts(info_dict, counts_iter, taxonomy, inc_anc=None,
                       rank=None, gene_map=None, ex_anc=None,
                       include_higher=True, cached=True, uid_used=None):
    """
    .. versionchanged:: 0.1.14
        added *cached* argument

    .. versionchanged:: 0.1.15
        added *uid_used* parameter

    .. versionchanged:: 0.2.0
        info_dict can be a function

    Reads sample counts, filtering and mapping them if requested. It's an
    example of the usage of the above functions.

    Arguments:
        info_dict (dict): dictionary that has *uid* as key and
            *(gene_id, taxon_id)* as value. In alternative a function that
            accepts a *uid* as sole argument and returns *(gene_id, taxon_id)*
        counts_iter (iterable): iterable that yields a *(uid, count)*
        taxonomy: taxonomy instance
        inc_anc (int, list): ancestor taxa to include
        rank (str): rank to which map the counts
        gene_map (dict): dictionary with the gene mappings
        ex_anc (int, list): ancestor taxa to exclude
        include_higher (bool): if False, any rank different than the requested
            one is discarded
        cached (bool): if *True*, the function will use
            :class:`mgkit.simple_cache.memoize` to cache some of the functions
            used
        uid_used (None, dict): an empty dictionary in which to store the *uid*
            that were assigned to each key of the returned pandas.Series. If
            *None*, no information is saved

    Returns:
        pandas.Series: array with MultiIndex *(gene_id, taxon_id)* with the
        filtered and mapped counts
    """
    if isinstance(info_dict, dict):
        if isinstance(info_dict[list(info_dict.keys())[0]], tuple):
            info_func = functools.partial(get_uid_info, info_dict)
        else:
            info_func = functools.partial(get_uid_info_ann, info_dict)
    else:
        info_func = info_dict

    tfilters = []

    if inc_anc is not None:
        if not isinstance(inc_anc, (list, set, tuple)):
            inc_anc = [inc_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=inc_anc,
                exclude=False,
                taxonomy=taxonomy
            )
        )
    if ex_anc is not None:
        if not isinstance(ex_anc, (list, set, tuple)):
            ex_anc = [ex_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=ex_anc,
                exclude=True,
                taxonomy=taxonomy
            )
        )

    if cached:
        tfilters = [
            mgkit.simple_cache.memoize(tfilter)
            for tfilter in tfilters
        ]

    gmapper = None
    tmapper = None

    if rank is not None:
        tmapper = functools.partial(
            map_taxon_id_to_rank,
            taxonomy,
            rank,
            include_higher=include_higher
        )
        if cached:
            tmapper = mgkit.simple_cache.memoize(tmapper)

    if gene_map is not None:
        gmapper = functools.partial(
            map_gene_id_to_map,
            gene_map
        )
        if cached:
            gmapper = mgkit.simple_cache.memoize(gmapper)

    series = map_counts(
        filter_counts(
            counts_iter,
            info_func,
            gfilters=None,
            tfilters=tfilters
        ),
        info_func,
        gmapper=gmapper,
        tmapper=tmapper,
        uid_used=uid_used
    )

    return series


def load_sample_counts_to_taxon(info_func, counts_iter, taxonomy, inc_anc=None,
                                rank=None, ex_anc=None, include_higher=True,
                                cached=True, uid_used=None):
    """
    .. versionadded:: 0.1.14

    .. versionchanged:: 0.1.15
        added *uid_used* parameter

    Reads sample counts, filtering and mapping them if requested. It's a
    variation of :func:`load_sample_counts`, with the counts being mapped only
    to each specific taxon. Another difference is the absence of any assumption
    on the first parameter. It is expected to return a (gene_id, taxon_id)
    tuple.

    Arguments:
        info_func (callable): any callable that accept an *uid* as the only
            parameter and and returns *(gene_id, taxon_id)* as value
        counts_iter (iterable): iterable that yields a *(uid, count)*
        taxonomy: taxonomy instance
        inc_anc (int, list): ancestor taxa to include
        rank (str): rank to which map the counts
        ex_anc (int, list): ancestor taxa to exclude
        include_higher (bool): if False, any rank different than the requested
            one is discarded
        cached (bool): if *True*, the function will use
            :class:`mgkit.simple_cache.memoize` to cache some of the functions
            used
        uid_used (None, dict): an empty dictionary in which to store the *uid*
            that were assigned to each key of the returned pandas.Series. If
            *None*, no information is saved

    Returns:
        pandas.Series: array with Index *taxon_id* with the filtered and mapped
        counts
    """
    tfilters = []

    if inc_anc is not None:
        if not isinstance(inc_anc, (list, set, tuple)):
            inc_anc = [inc_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=inc_anc,
                exclude=False,
                taxonomy=taxonomy
            )
        )
    if ex_anc is not None:
        if not isinstance(ex_anc, (list, set, tuple)):
            ex_anc = [ex_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=ex_anc,
                exclude=True,
                taxonomy=taxonomy
            )
        )

    if cached:
        tfilters = [
            mgkit.simple_cache.memoize(tfilter)
            for tfilter in tfilters
        ]

    tmapper = None

    if rank is not None:
        tmapper = functools.partial(
            map_taxon_id_to_rank,
            taxonomy,
            rank,
            include_higher=include_higher
        )
        if cached:
            tmapper = mgkit.simple_cache.memoize(tmapper)

    series = map_counts(
        filter_counts(
            counts_iter,
            info_func,
            gfilters=None,
            tfilters=tfilters
        ),
        info_func,
        tmapper=tmapper,
        index='taxon',
        uid_used=uid_used
    )

    return series


def load_sample_counts_to_genes(info_func, counts_iter, taxonomy, inc_anc=None,
                                gene_map=None, ex_anc=None, cached=True,
                                uid_used=None):
    """
    .. versionadded:: 0.1.14

    .. versionchanged:: 0.1.15
        added *uid_used* parameter

    Reads sample counts, filtering and mapping them if requested. It's a
    variation of :func:`load_sample_counts`, with the counts being mapped only
    to each specific gene_id. Another difference is the absence of any
    assumption on the first parameter. It is expected to return a
    (gene_id, taxon_id) tuple.

    Arguments:
        info_func (callable): any callable that accept an *uid* as the only
            parameter and and returns *(gene_id, taxon_id)* as value
        counts_iter (iterable): iterable that yields a *(uid, count)*
        taxonomy: taxonomy instance
        inc_anc (int, list): ancestor taxa to include
        rank (str): rank to which map the counts
        gene_map (dict): dictionary with the gene mappings
        ex_anc (int, list): ancestor taxa to exclude
        cached (bool): if *True*, the function will use
            :class:`mgkit.simple_cache.memoize` to cache some of the functions
            used
        uid_used (None, dict): an empty dictionary in which to store the *uid*
            that were assigned to each key of the returned pandas.Series. If
            *None*, no information is saved

    Returns:
        pandas.Series: array with Index *gene_id* with the filtered and mapped
        counts
    """
    tfilters = []

    if inc_anc is not None:
        if not isinstance(inc_anc, (list, set, tuple)):
            inc_anc = [inc_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=inc_anc,
                exclude=False,
                taxonomy=taxonomy
            )
        )
    if ex_anc is not None:
        if not isinstance(ex_anc, (list, set, tuple)):
            ex_anc = [ex_anc]
        tfilters.append(
            functools.partial(
                tx_filters.filter_by_ancestor,
                filter_list=ex_anc,
                exclude=True,
                taxonomy=taxonomy
            )
        )

    if cached:
        tfilters = [
            mgkit.simple_cache.memoize(tfilter)
            for tfilter in tfilters
        ]

    if gene_map is not None:
        gmapper = functools.partial(
            map_gene_id_to_map,
            gene_map
        )
        if cached:
            gmapper = mgkit.simple_cache.memoize(gmapper)

    series = map_counts(
        filter_counts(
            counts_iter,
            info_func,
            gfilters=None,
            tfilters=tfilters
        ),
        info_func,
        gmapper=gmapper,
        index='gene',
        uid_used=uid_used
    )

    return series


def load_deseq2_results(file_name, taxon_id=None):
    """
    .. versionadded:: 0.1.14

    Reads a CSV file output with DESeq2 results, adding a taxon_id to the index
    for concatenating multiple results from different taxonomic groups.

    Arguments:
        file_name (str): file name of the CSV
    Returns:
        pandas.DataFrame: a MultiIndex DataFrame with the results

    """
    dataframe = pandas.DataFrame.from_csv(file_name)

    dataframe = dataframe.rename(
        index=dict(
            (gene_id, (gene_id, taxon_id))
            for gene_id in dataframe.index
        )
    )
    dataframe.index.names = ['gene_id', 'taxon_id']

    return dataframe


def map_counts_to_category(counts, gene_map, nomap=False, nomap_id='NOMAP'):
    """
    Used to map the counts from a certain gene identifier to another. Genes
    with no mappings are not counted, unless *nomap=True*, in which case they
    are counted as *nomap_id*.

    Arguments:
        counts (iterator): an iterator that yield a tuple, with the first value
            being the gene_id and the second value the count for it
        gene_map (dictionary): a dictionary whose keys are the gene_id yield by
            *counts* and the values are iterable of mapping identifiers
        nomap (bool): if False, counts for genes with no mappings in *gene_map*
            are discarded, if True, they a counted as *nomap_id*
        nomap_id (str): name of the mapping for genes with no mappings

    Returns:
        pandas.Series: mapped counts
    """

    newcounts = {}
    for gene_id, count in counts:
        try:
            map_ids = gene_map[gene_id]
        except KeyError:
            continue

        if nomap and (not map_ids):
            map_ids = [nomap_id]

        for map_id in map_ids:
            try:
                newcounts[map_id] += count
            except KeyError:
                newcounts[map_id] = count

    return pandas.Series(newcounts)


def load_counts_from_gff(annotations, elem_func=lambda x: x.uid, sample_func=None, nozero=True):
    """
    .. versionadded:: 0.2.5

    Loads counts for each annotations that are stored into the annotation
    *counts_* attributes. Annotations with a total of 0 counts are skipped by
    default (nozero=True), the row index is set to the *uid* of the annotation
    and the column to the sample name. The functions used to transform the
    indices expect the annotation (for the row, *elem_func*) and the sample
    name (for the column, *sample_func*).

    Arguments:
        annotations (iter): iterable of annotations
        elem_func (func): function that accepts an annotation and return a
            str/int for a Index or a tuple for a MultiIndex, defaults to
            returning the *uid* of the annotation
        sample_func (func, None): function that accepts the sample name and
            returns tuple for a MultiIndex. Defaults to *None* so no
            transformation is performed
        nozero (bool): if *True*, annotations with no counts are skipped
    """

    count_dict = {}

    for annotation in annotations:
        counts = pandas.Series(annotation.counts)

        if nozero and (counts.sum() == 0):
            continue

        count_dict[elem_func(annotation)] = counts

    # LOG.debug('Returning DataFrame')

    count_dict = pandas.DataFrame.from_dict(count_dict, orient='index')
    if sample_func is not None:
        count_dict.columns = pandas.MultiIndex.from_tuples(
            [
                sample_func(column)
                for column in count_dict.columns
            ]
        )

    return count_dict


def from_gff(annotations, samples, ann_func=None, sample_func=None):
    """
    .. versionadded:: 0.3.1

    Loads count data from a GFF file, only for the requested samples. By
    default the function returns a DataFrame where the index is the *uid* of
    each annotation and the columns the requested samples.

    This can be customised by supplying *ann_func* and *sample_func*.
    *sample_func* is a function that accept a sample name and is expected to
    return a string or a tuple. This will be used to change the columns in the
    DataFrame. *ann_func* must accept an :class:`mgkit.io.gff.Annotation`
    instance and return an iterable, with each iteration yielding either a
    single element or a tuple (for a MultiIndex DataFrame), each element
    yielded will have the count of that annotation added to.

    Arguments:
        annotation (iterable): iterable yielding annotations
        samples (iterable): list of samples to keep
        ann_func (func): function used to customise the output
        sample_func (func): function to customise the column elements

    Returns:
        DataFrame: dataframe with the count data, columns are the samples and
        rows the annotation counts (unless mapped with *ann_func*)

    Exmples:
        Assuming we have a list of *annotations* and sample SAMPLE1 and SAMPLE2
        we can obtain the count table for all annotations with this

        >>> from_gff(annotations, ['SAMPLE1', 'SAMPLE2'])

        Assuming we want to group the samples, for example treatment1,
        treatment2 and control1, control2 into a MultiIndex DataFrame column

        >>> sample_func = lambda x: ('T' if x.startswith('t') else 'C', x)
        >>> from_gff(annotations, ['treatment1', 'treatment2', 'control1', 'control2'], sample_func=sample_func)

        Annotations can be mapped to other levels for example instead of using
        the *uid* that is the default, it can be mapped to the gene_id, taxon_id
        information that is included in the annotation, resulting in a MultiIndex
        index for the rows, with (gene_id, taxon_id) as key.

        >>> ann_func = lambda x: [(x.gene_id, x.taxon_id)]
        >>> from_gff(annotations, ['SAMPLE1', 'SAMPLE2'], ann_func=ann_func)
    """
    if ann_func is None:
        ann_func = lambda x: (x.uid, )

    if sample_func is None:
        sample_func = lambda x: x

    counters = {
        sample_func(sample): Counter()
        for sample in samples
    }

    for annotation in annotations:
        for sample, count in viewitems(annotation.counts):
            sample = sample_func(sample)
            if sample not in counters:
                continue
            for ann_id in ann_func(annotation):
                counters[sample][ann_id] += count

    return pandas.DataFrame(counters)
