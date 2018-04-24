import pytest

from mgkit.utils.common import average_length, between, union_range, \
    range_intersect, union_ranges, complement_ranges, ranges_length, \
    apply_func_window


@pytest.mark.parametrize(
    "coords,result",
    [
        ((1, 200, 3, 132), 165.0),
        ((4, 2245, 3, 13223), 7731.5),
    ]
)
def test_avg_len(coords, result):
    assert average_length(*coords) == result


@pytest.mark.parametrize(
    "coords,result",
    (
        ((1, 0, 10), True),
        ((1, 2, 10), False),
        ((1, 1, 10), True),
        ((11, 1, 10), False),
        ((10, 1, 10), True),
    )
)
def test_between(coords, result):
    assert between(*coords) == result


@pytest.mark.parametrize(
    "coords,result",
    [
        ((1, 10, 2, 13), (1, 13)),
        ((1, 10, 10, 13), (1, 13)),
        ((10, 13, 1, 10), (1, 13)),
        ((1, 10, 11, 13), (1, 13)),
        ((10.0, 13.0, 1, 10), (1, 13)),
    ]
)
def test_union_range(coords, result):
    assert union_range(*coords) == result


@pytest.mark.parametrize(
    "coords,result",
    [
        ((10, 20, 19, 30), (19, 20)),
        ((10, 20, 15, 30), (15, 20)),
        ((10, 20, 10, 20), (10, 20)),
        ((10, 20, 12, 18), (12, 18)),
        # The next one is to test the result for no intersection
        ((10, 20, 30, 40), None)
    ]
)
def test_range_intersect(coords, result):
    assert range_intersect(*coords) == result


@pytest.mark.parametrize(
    "intervals,result",
    [
        ([(1, 2), (3, 7), (6, 12), (9, 17), (18, 20)], [(1, 20)]),
        ([(3, 7), (6, 12), (9, 17), (18, 20), (1, 2)], [(1, 20)]),
        ([(1, 2), (3, 7), (6, 12), (9, 14), (18, 20)], [(1, 14), (18, 20)])
    ]
)
def test_union_ranges(intervals, result):
    assert union_ranges(intervals) == result


@pytest.mark.parametrize(
    "intervals,end,result",
    [
        ([(1, 10), (11, 20), (25, 30)], 100, [(21, 24), (31, 100)]),
        ([(1, 10), (11, 20), (25, 30)], None, [(21, 24)]),
        ([(0, 2), (3, 17), (18, 20)], None, []),
        ([(0, 2), (3, 17), (18, 20)], 100, [(21, 100)])
    ]
)
def test_complement_ranges(intervals, end, result):
    assert complement_ranges(intervals, end) == result


@pytest.mark.parametrize(
    "coords,result",
    [
        ([(1, 10), (12, 25), (25, 30)], 30),
        # They overlap by one,, so the sum will be 26, instead of 25
        # Tested to make sure we can rely on this behaviour
        # as per warning in the documentation
        ([(1, 10), (10, 25)], 26),
    ]
)
def test_range_intersect(coords, result):
    assert ranges_length(coords) == result


@pytest.mark.parametrize(
    "func,data,window,step,result",
    [
        [sum, range(10), 10, 0, [45]],
        [sum, range(10), 10, 9, [45, 9]],
        [sum, range(10), 5, 5, [10, 35]],
    ]
)
def test_apply_func_window(func, data, window, step, result):
    reduction = list(
        apply_func_window(func, data, window, step)
    )

    assert reduction == result
