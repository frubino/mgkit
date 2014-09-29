"""
.. versionadded:: 0.1.13

Misc function for count data
"""

import logging
import itertools
import pandas
import functools
from mgkit.filter import taxon as tx_filters
from mgkit.io import open_file

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


def load_htseq_counts(f_handle):
    """
    Loads an HTSeq-count result file

    Arguments:

        f_handle (file or str): file handle or string with file name

    Yields:
        tuple: first element is the gene_id and the second is the count

    """

    LOG.info("Loading HTSeq-count file %s", str(f_handle))

    for line in open_file(f_handle, 'r'):
        gene_id, count = line.rstrip().split('\t')

        if line.startswith('__') or (gene_id in SKIP):
            continue

        yield gene_id, int(count)


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

    iterator = itertools.izip_longest(
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


def map_counts(counts_iter, info_func, gmapper=None, tmapper=None):
    """
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
            key = (map_id, taxon_id)

            try:
                mapped_counts[key] += count
            except KeyError:
                mapped_counts[key] = count

    return pandas.Series(mapped_counts)


def map_taxon_id_to_rank(taxonomy, rank, taxon_id, include_higher=True):
    """
    Maps a *taxon_id* to the request taxon rank. Returns *None* if
    *include_higher* is False and the found rank is not the one requested.

    Internally uses :meth:`mgkit.taxon.UniprotTaxonomy.get_ranked_taxon`

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
                       include_higher=True):
    """
    Reads sample counts, filtering and mapping them if requested. It's an
    example of the usage of the above functions.

    Arguments:
        info_dict (dict): dictionary that has *uid* as key and
            *(gene_id, taxon_id)* as value
        counts_iter (iterable): iterable that yields a *(uid, count)*
        taxonomy: taxonomy instance
        inc_anc (int, list): ancestor taxa to include
        rank (str): rank to which map the counts
        gene_map (dict): dictionary with the gene mappings
        ex_anc (int, list): ancestor taxa to exclude
        include_higher (bool): if False, any rank different than the requested
            one is discarded

    Returns:
        pandas.Series: array with MultiIndex *(gene_id, taxon_id)* with the
        filtered and mapped counts
    """
    if isinstance(info_dict[info_dict.keys()[0]], tuple):
        info_func = functools.partial(get_uid_info, info_dict)
    else:
        info_func = functools.partial(get_uid_info_ann, info_dict)

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

    gmapper = None
    tmapper = None

    if rank is not None:
        tmapper = functools.partial(
            map_taxon_id_to_rank,
            taxonomy,
            rank,
            include_higher=include_higher
        )

    if gene_map is not None:
        gmapper = functools.partial(
            map_gene_id_to_map,
            gene_map
        )

    series = map_counts(
        filter_counts(
            counts_iter,
            info_func,
            gfilters=None,
            tfilters=tfilters
        ),
        info_func,
        gmapper=gmapper,
        tmapper=tmapper
    )

    return series
