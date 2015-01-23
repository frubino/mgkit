"""
.. versionadded:: 0.1.14

Contains code related to colour
"""
import logging


LOG = logging.getLogger(__name__)


def float_to_hex_color(r, g, b):
    """
    .. versionadded:: 0.1.14

    Converts RGB float values to Hexadecimal value string
    """
    convert = lambda x: int(x * 255)

    return "#{0:02x}{1:02x}{2:02x}".format(convert(r), convert(g), convert(b))

__all__ = ['float_to_hex_color']
