from nose import SkipTest


def skip_test(func, msg=''):
    def wrapper(*args, **kwd):
        raise SkipTest(msg)

    # This makes sure that the test is run, otherwise the name of the
    # decorated function will be "skip_test", which is in the excluded
    # list
    wrapper.__name__ = func.__name__

    return wrapper
