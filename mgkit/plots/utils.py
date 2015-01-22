"""
.. versionadded:: 0.1.14

Misc code

"""
import logging

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

LOG = logging.getLogger(__name__)


def get_single_figure(dpi=300, figsize=(10, 20), aspect='auto'):
    """
    .. versionchanged:: 0.1.14
        added *aspect* parameter

    Simple wrapper to init a single figure

    Arguments:
        dpi (int): dpi used for the figure
        figsize (tuple): size of the figure in inches
        aspect (str, float): aspect ratio to be passed to figure.add_subplot

    Returns:
        tuple: the figure and axes objects
    """
    fig = plt.figure(dpi=dpi, figsize=figsize)
    ax = fig.add_subplot(111, aspect=aspect)
    return fig, ax


def get_grid_figure(rows, cols, dpi=300, figsize=(10, 20), **kwd):
    """
    .. versionadded:: 0.1.13

    Simple wrapper to init a GridSpec figure

    Arguments:
        rows (int): number of rows
        columns (int): number of columns
        dpi (int): dpi used for the figure
        figsize (tuple): size of the figure in inches

    Returns:
        tuple: the figure and axes objects
    """
    fig = plt.figure(dpi=dpi, figsize=figsize)
    gs = GridSpec(rows, cols, **kwd)
    return fig, gs


__all__ = ['get_grid_figure', 'get_single_figure']
