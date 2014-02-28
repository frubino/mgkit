from nose.tools import *

from mgkit.utils.common import average_length, between


def test_avg_len1():
    eq_(average_length(1, 200, 3, 132), 165.0)


def test_avg_len2():
    eq_(average_length(4, 2245, 3, 13223), 7731.5)


def test_between1():
    eq_(between(1, 0, 10), True)


def test_between2():
    eq_(between(1, 2, 10), False)


def test_between3():
    eq_(between(1, 1, 10), True)


def test_between4():
    eq_(between(11, 1, 10), False)


def test_between5():
    eq_(between(10, 1, 10), True)
