"""
Access some R functions. requires :mod:`rpy2`
"""

#to get python 3.x behavior on relative imports - it's absolute unless specified
#by a dot - otherwise current dir module is loaded
from __future__ import absolute_import
import rpy2.robjects as robjects
import pandas


def correct_pvalues(pvalues, method='BH'):
    """
    Correct pvalues using R.

    :param pandas.Series pvalues: :class:`pandas.Series` instance containing the
        pvalues
    :param str method: string passed to the R function to specify the correction
        method

    :return: a :class:`pandas.Series` with the same index as pvalues
    """

    padj = robjects.r('p.adjust')

    corr_pvalues = padj(robjects.FloatVector(pvalues), method)
    corr_pvalues = pandas.Series(corr_pvalues, index=pvalues.index)

    return corr_pvalues
