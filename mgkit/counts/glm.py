"""
.. versionadded:: 0.3.3

GLM models with metagenomes and metatranscriptomes. Experimental
"""

try:
    import statsmodels.api as sm
    from scipy import interpolate
    from scipy import optimize
    import pandas as pd
except ImportError:
    raise DependencyError('statsmodels, scipy, pandas')


def lowess_ci_bootstrap(endog, exog, num=100, frac=.2, it=3, alpha=.05):
    data = pd.DataFrame(
        {
            'endog': endog,
            'exog': exog
        }
    ).sort_values(by=['exog', 'endog'], ascending=True)

    boots = []
    for idx in xrange(num):
        sample = data.sample(n=data.index.size, replace=True)
        lw = sm.nonparametric.lowess(
            sample.endog,
            sample.exog,
            is_sorted=False,
            frac=frac,
            it=it,
            return_sorted=True
        )
        boots.append(pd.DataFrame({
            'endog': lw[:, 1],
            'exog': lw[:, 0]
        }))

    boots = pd.concat(boots)
    q1 = boots.groupby('exog').quantile(
        alpha, interpolation='nearest'
    ).sort_index()
    q2 = boots.groupby('exog').quantile(
        1 - alpha, interpolation='nearest'
    ).sort_index()
    return q1, q2


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
    )
