"""
.. versionadded:: 0.3.3

GLM models with metagenomes and metatranscriptomes. Experimental
"""

from __future__ import division
from builtins import range
import statsmodels.api as sm
from scipy import interpolate, optimize
import pandas as pd


def lowess_ci_bootstrap(endog, exog, num=100, frac=.2, it=3, alpha=.05,
                        delta=0., min_value=10**-3, kind='slinear'):
    """
    .. versionchanged: 0.4.0
        the interpolation of the confidence intervals uses the passed `kind`

    Bootstraps a lowess for the dependent (`endog`) and indipendent (`exog`)
    arguments.

    Arguments:
        endog (array): indipendent variable (Y)
        exog (array): indipendent variable (X)
        num (int): number of iterations for the bootstrap
        frac (float): fraction of the array to use when fitting
        it (int): number of iterations used to fit the lowess
        alpha (float): confidence intervals for the bootstrap
        delta (float): passed to :func:`statsmodels.api.nonparametric.lowess`
        min_value (float): minimum value for the function to avoid out of
            bounds
        kind (str): type of interpolation passed to
            :func:`scipy.interpolate.interp1d`

    Returns:
        tuple: the first element is the function describing the lowest
        confidence interval, the second element is for the highest confidence
        interval and the last one for the mean

    .. note::

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
    q1 = q1.endog[q1.endog > 0].fillna(min_value)
    q1 = interpolate.interp1d(q1.index, q1, kind=kind,
                              fill_value='extrapolate')

    q2 = boots.groupby('exog').quantile(
        1 - alpha, interpolation='nearest'
    ).sort_index()
    q2 = q2.endog[q2.endog > 0].fillna(min_value)
    q2 = interpolate.interp1d(q2.index, q2, kind=kind,
                              fill_value='extrapolate')

    m = fit_lowess_interpolate(data.endog, data.exog, frac=frac, it=it,
                               kind=kind)

    return q1, q2, m


def fit_lowess_interpolate(endog, exog, frac=.2, it=3, kind='slinear'):
    """
    Fits a lowess for the passed `endog` (Y) and `exog` (X) and returns an
    interpolated function that describes it. The first 4 arguments are passed
    to :func:`statsmodels.api.sm.nonparametric.lowess`, while the last one is
    passed to :func:`scipy.interpolate.interp1d`

    Arguments:
        endog (array): array of the dependent variable (Y)
        exog (array): array of the indipendent variable (X)
        frac (float): fraction of the number of elements to use when fitting
            (0.0-1.0)
        it (int): number of iterations to fit the lowess
        kind (str): type of interpolation to use

    Returns:
        func: interpolated function representing the lowess fitted from the
        data passed
    """

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

    Arguments:
        mu (float): mean to calculate the alphas for
        func (func): function that returns the variace of the mean
        min_alpha (float): value of alpha if the `func` goes out of bounds

    Returns:
        float: value of alpha for the passed mean
    """

    alpha = (func(mu) - mu) / (mu ** 2)

    return min_alpha if alpha <= 0 else alpha


def optimise_alpha_scipy_function(args, formula, data, criterion='aic'):
    """
    .. versionadded:: 0.4.0


    """
    alpha, = args
    model = sm.formula.glm(
        formula,
        data=data,
        family=sm.families.NegativeBinomial(alpha=alpha)
    ).fit()
    return getattr(model, criterion)


def optimise_alpha_scipy(formula, data, mean_func, q1_func, q2_func):
    """
    .. versionadded:: 0.4.0

    Used to find an optimal *alpha* parameter for the Negative Binomial
    distribution used in `statsmodels`, using the lowess functions from
    :func:`lowess_ci_bootstrap`.

    Arguments:
        formula (str): the formula used for the regression
        data (DataFrame): DataFrame for regression
        mean_func (func): function for the mean :func:`lowess_ci_bootstrap`
        q1_func (func): function for the q1 :func:`lowess_ci_bootstrap`
        q2_func (func): function for the q2 :func:`lowess_ci_bootstrap`

    Returns:
        float: *alpha* value for the Negative Binomial
    """
    data_mean = data.mean()[0]
    return optimize.minimize(
        optimise_alpha_scipy_function,
        [variance_to_alpha(data_mean, mean_func)],
        args=(formula, data),
        bounds=[(variance_to_alpha(data_mean, q1_func), variance_to_alpha(data_mean, q2_func))]
    ).x[0]
