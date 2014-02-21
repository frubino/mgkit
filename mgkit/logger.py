"Module configuring log information"

import logging
import sys


def config_log(level=logging.DEBUG, output=sys.stderr):
    """
    Minimal configuration of :mod`logging` module, default to debug level and
    the output is printed to standard error

    :param int level: logging level
    :param file output: file to which write the log
    """
    log_handler = logging.StreamHandler(output)
    if level == logging.DEBUG:
        fmt = logging.Formatter(
            fmt="%(asctime)s - %(levelname) 7s - %(name)s->%(funcName)s: " +
                "%(message)s")
    else:
        fmt = logging.Formatter(
            fmt="%(levelname)s - %(name)s: %(message)s")
    log_handler.setFormatter(fmt)
    logging.getLogger().addHandler(log_handler)
    logging.getLogger().setLevel(level)
