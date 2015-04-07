"""
Utility functions for workflows
"""
from __future__ import print_function
import logging
import argparse
import sys
import mgkit

LOG = logging.getLogger(__name__)


class CiteAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help='Show citation for the framework'):
        super(CiteAction, self).__init__(
            option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        mgkit.cite(file=sys.stderr)
        setattr(namespace, self.dest, values)
        parser.exit(0)


def add_basic_options(parser):
    """
    Adds verbose and version options to the option parser
    """
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_const',
        const=logging.DEBUG,
        default=logging.INFO,
        help='more verbose'
    )
    parser.add_argument(
        '--cite',
        action=CiteAction
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {0}'.format(mgkit.__VERSION__)
    )


def exit_script(message, ret_value):
    LOG.critical(message)
    sys.exit(ret_value)
