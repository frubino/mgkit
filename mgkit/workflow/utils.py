"""
Utility functions for workflows
"""
from __future__ import print_function
import logging
import argparse
import sys
import click
import mgkit

LOG = logging.getLogger(__name__)


def cite_callback(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(mgkit.CITE)
    ctx.exit()


cite_option = click.option('--cite', is_flag=True, callback=cite_callback,
                           expose_value=False, is_eager=True)


class CiteAction(argparse.Action):
    """
    Argparse action to print the citation, using the :func:`mgkit.cite`
    function.
    """
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
        mgkit.cite(file_handle=sys.stderr)
        setattr(namespace, self.dest, values)
        parser.exit(0)


class PrintManAction(argparse.Action):
    """
    .. versionadded:: 0.2.6

    Argparse action to print the manual
    """
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help='Show the script manual',
                 manual=''):
        super(PrintManAction, self).__init__(
            option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)
        self.manual = manual

    def __call__(self, parser, namespace, values, option_string=None):
        print(self.manual, file=sys.stdout)
        setattr(namespace, self.dest, values)
        parser.exit(0)


def add_basic_options(parser, manual=''):
    """
    .. versionchanged:: 0.2.6
        added *quiet* option

    Adds verbose and version options to the option parser
    """
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-v',
        '--verbose',
        action='store_const',
        const=logging.DEBUG,
        default=logging.INFO,
        help='more verbose - includes debug messages',
        dest='verbose'
    )
    group.add_argument(
        '--quiet',
        action='store_const',
        const=logging.ERROR,
        help='less verbose - only error and critical messages',
        dest='verbose'
    )
    parser.add_argument(
        '--cite',
        action=CiteAction
    )
    parser.add_argument(
        '--manual',
        action=PrintManAction,
        manual=manual
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {0}'.format(mgkit.__VERSION__)
    )


def exit_script(message, ret_value):
    """
    Used to exit the script with a return value
    """
    LOG.critical(message)
    sys.exit(ret_value)
