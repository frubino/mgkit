"Module configuring log information"

import logging
import sys
import os.path

DEBUG_FMT = "%(asctime)s - %(levelname) 7s - %(name)s->%(funcName)s: %(message)s"
INFO_FMT = "%(levelname)s - %(name)s: %(message)s"


def config_log(level=logging.DEBUG, output=sys.stderr):
    """
    Minimal configuration of :mod`logging` module, default to debug level and
    the output is printed to standard error

    :param int level: logging level
    :param file output: file to which write the log
    """
    log_handler = logging.StreamHandler(output)

    if level == logging.DEBUG:
        fmt = logging.Formatter(fmt=DEBUG_FMT)
    else:
        fmt = logging.Formatter(fmt=INFO_FMT)

    log_handler.setFormatter(fmt)
    logging.getLogger().addHandler(log_handler)
    logging.getLogger().setLevel(level)


def config_log_to_file(level=logging.DEBUG, output=None):
    """
    .. versionadded:: 0.1.14

    Minimal configuration of :mod`logging` module, default to debug level and
    the output is printed to script name, using `sys.argv[0]`.

    :param int level: logging level
    :param file output: file to which write the log
    """

    if output is None:
        output = os.path.join(
            os.path.splitext(sys.argv[0])[0],
            '.log'
        )

    log_handler = logging.FileHandler(filename=output, mode='w', delay=True)

    if level == logging.DEBUG:
        fmt = logging.Formatter(fmt=DEBUG_FMT)
    else:
        fmt = logging.Formatter(fmt=INFO_FMT)

    log_handler.setFormatter(fmt)
    logging.getLogger().addHandler(log_handler)
    logging.getLogger().setLevel(level)
