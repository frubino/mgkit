"""
Scaling functions for counts
"""

import scipy
import numpy
import pandas


def scale_factor_deseq(dataframe):
    """
    Returns the scale factor according to he deseq paper. The columns of the
    dataframe are the samples.

    size factor :math:`\\hat{s}_{j}` from deseq paper

    .. math::

        \\hat{s}_{j} =
        \\frac
            {k_{ij}}
            {
                \\left (
                \prod_{v=1}^{m}
                    k_{iv}
                \\right )^{1/m}
           }

    """
    #calc the genes geometric mean over all samples
    gmean = dataframe.apply(scipy.stats.gmean, axis=1)
    #keeps only the genes whose geometric mean is > 0
    gmean = gmean[gmean > 0]

    sample_factors = {}

    #calc the scaling factor for each sample
    for sample, genes in dataframe.iterkv():

        scale_factor = numpy.median(genes.loc[gmean.index] / gmean)

        sample_factors[sample] = scale_factor

    return pandas.Series(sample_factors)


def scale_deseq(dataframe):
    """
    Scale a dataframe using the deseq scaling. Uses :func:`scale_factor_deseq`
    """
    scale_factors = scale_factor_deseq(dataframe)

    return dataframe / scale_factors
