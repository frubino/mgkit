from nose.tools import *
from nose import SkipTest

try:
    from mgkit.utils import r_func
except ImportError:
    raise SkipTest('Cannot import r_func, rpy2 not installed?')
import pandas
import numpy


def setup_random_vector():
    return pandas.Series(numpy.random.random_sample(1000))


def test_correct_pvalues():
    pvalues = setup_random_vector()
    assert len(r_func.correct_pvalues(pvalues)) == len(pvalues)
