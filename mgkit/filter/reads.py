"""
Some test functions to filter sequences
"""
from __future__ import division

import numpy
# from numpy import sum
# import pandas
# from scipy.stats import itemfreq
# import numexpr as ne

# ne.set_num_threads(8)


# #direct use: fourth (1.6ms loop)
# #in trim_by_ee: fourth
# def expected_error_rate_1(qualities):
#     qualities = pandas.DataFrame(itemfreq(qualities))

#     return qualities.prod(axis=1).sum()


#direct use: first (0.033ms loop)
#in trim_by_ee: first (15.9s)
def expected_error_rate(qualities):
    """
    Calculate the expected error rate for an array of qualities (converted to
    probabilities).
    """
    unique, inv = numpy.unique(qualities, return_inverse=True)
    return sum(unique * numpy.bincount(inv))


# #direct use: third (0.463ms loop)
# #in trim_by_ee: third (65.7s)
# def expected_error_rate_3(qualities):
#     return sum(
#         sum(1 for value in qualities if value == unique) * unique
#         for unique in numpy.unique(qualities)
#     )


# #direct use: second (0.109ms loop)
# #in trim_by_ee: second (27.7s)
# def expected_error_rate_4(qualities):
#     return sum(
#         len(qualities[qualities == unique]) * unique
#         for unique in numpy.unique(qualities)
#     )


def trim_by_ee(qualities, min_length=50, threshold=0.5, chars=True, base=33):
    """
    Trim a sequence based on the expected error rate.
    """
    if chars:
        qualities = numpy.fromiter(
            (ord(score) for score in qualities),
            dtype=int
        ) - base
    # qualities = ne.evaluate("10 ** ((-qualities) / 10)")
    qualities = 10 ** ((-qualities) / 10)
    # print qualities
    for idx in xrange(min_length, qualities.size):
        # print expected_error_rate(qualities[:idx])
        if expected_error_rate(qualities[:idx]) > threshold:
            return idx
    return qualities.size
