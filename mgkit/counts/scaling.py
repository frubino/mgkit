"""
Scaling functions for counts
"""

from scipy import stats
import numpy
import pandas


def scale_factor_deseq(dataframe):
    """
    .. versionadded:: 0.1.13

    Returns the scale factor according to he deseq paper. The columns of the
    dataframe are the samples.

    size factor :math:`\\hat{s}_{j}` for sample *j* (from DESeq paper).

    .. math::

        \\hat{s}_{j} = median_{i} (
        \\frac
            {k_{ij}}
            {
                \\left (
                \\prod_{v=1}^{m}
                    k_{iv}
                \\right )^{1/m}
           }
        )

    """
    # calc the genes geometric mean over all samples
    gmean = dataframe.apply(stats.gmean, axis=1)
    # keeps only the genes whose geometric mean is > 0
    gmean = gmean[gmean > 0]

    sample_factors = {}

    # calc the scaling factor for each sample
    for sample, genes in dataframe.iteritems():

        scale_factor = numpy.median(genes.loc[gmean.index] / gmean)

        sample_factors[sample] = scale_factor

    return pandas.Series(sample_factors)


def scale_deseq(dataframe):
    """
    .. versionadded:: 0.1.13

    Scale a dataframe using the deseq scaling. Uses :func:`scale_factor_deseq`
    """
    scale_factors = scale_factor_deseq(dataframe)

    return dataframe / scale_factors


def scale_rpkm(dataframe, gene_len):
    """
    .. versionadded:: 0.1.14

    Perform an RPKM scaling of the pandas dataframe/series supplied using the
    *gene_len* series containing the gene sizes for all elements of *dataframe*

    .. math::

        RPKM =\\frac {10^{9} \\cdot C} {N \\cdot L}

    """

    gene_len = gene_len[dataframe.index]
    tot_reads = dataframe.sum().sum()

    return (10 ** 9) * dataframe.div(gene_len * tot_reads, axis='index')
