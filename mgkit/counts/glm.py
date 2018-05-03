"""
.. versionadded:: 0.3.3

GLM models with metagenomes and metatranscriptomes. Experimental
"""

from __future__ import division
from builtins import range
import statsmodels.api as sm
from scipy import interpolate
from scipy import optimize
import pandas as pd


def lowess_ci_bootstrap(endog, exog, num=100, frac=.2, it=3, alpha=.05,
                        delta=0., min_alpha=10**-3, kind='slinear'):
    """
    Performance increase with the value of *delta*.
    """
    data = pd.DataFrame(
        {
            'endog': endog,
            'exog': exog
        }
    ).sort_values(by=['exog', 'endog'], ascending=True)

    boots = []
    for idx in range(num):
        sample = data.sample(n=data.index.size, replace=True)
        lw = sm.nonparametric.lowess(
            sample.endog,
            sample.exog,
            is_sorted=False,
            frac=frac,
            it=it,
            return_sorted=True,
            delta=delta
        )
        boots.append(pd.DataFrame({
            'endog': lw[:, 1],
            'exog': lw[:, 0]
        }))

    alpha = alpha / 2
    boots = pd.concat(boots)
    q1 = boots.groupby('exog').quantile(
        alpha, interpolation='nearest'
    ).sort_index()
    q1 = q1.endog[q1.endog > 0].fillna(min_alpha)
    q1 = interpolate.interp1d(q1.index, q1, kind='slinear', fill_value='extrapolate')

    q2 = boots.groupby('exog').quantile(
        1 - alpha, interpolation='nearest'
    ).sort_index()
    q2 = q2.endog[q2.endog > 0].fillna(min_alpha)
    q2 = interpolate.interp1d(q2.index, q2, kind='slinear', fill_value='extrapolate')

    m = fit_lowess_interpolate(data.endog, data.exog, frac=frac, it=it, kind=kind)

    return q1, q2, m

def fit_lowess_interpolate(endog, exog, frac=.2, it=3, kind='slinear'):
    data = pd.DataFrame(
        {
            'endog': endog,
            'exog': exog
        }
    ).sort_values(by=['exog', 'endog'], ascending=True)

    lw = sm.nonparametric.lowess(
        data.endog,
        data.exog,
        frac=frac,
        it=it,
    )

    return interpolate.interp1d(
        lw[:, 0],
        lw[:, 1],
        kind=kind,
        fill_value='extrapolate'
    )


def variance_to_alpha(mu, func, min_alpha=10**-3):
    """
    Based on the variance defined in the Negative Binomial in statsmodels

    var = mu + alpha * (mu ** 2)
    """

    alpha = (func(mu) - mu) / (mu ** 2)

    return min_alpha if alpha <= 0 else alpha
