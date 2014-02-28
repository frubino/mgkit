"""
Utility functions for workflows
"""
import logging
import mgkit


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
        '--version',
        action='version',
        version='%(prog)s {0}'.format(mgkit.__VERSION__)
    )
