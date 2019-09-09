import pytest
import pandas
import numpy
from mgkit.counts.scaling import scale_factor_deseq, scale_deseq


def test_scale_scale_factor_deseq():
    array = numpy.random.negative_binomial(100, .01, size=1000)
    dataframe = pandas.DataFrame(
        {1: array, 2: array * 1.5, 3: array * 2}
    )
    sf = scale_factor_deseq(dataframe)
    assert sf[1] < sf[2] < sf[3]


def test_scale_scale_deseq():
    array = numpy.random.negative_binomial(100, .01, size=1000)
    dataframe = pandas.DataFrame(
        {1: array, 2: array * 1.5, 3: array * 2}
    )
    sdf = scale_deseq(dataframe)
    sf = scale_factor_deseq(dataframe)
    assert (sdf * sf).round(0).iloc[0, :].sum() == dataframe.round(0).iloc[0, :].sum()
