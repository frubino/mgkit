"""
.. versionadded:: 0.1.15

Module to plot relative abundances in a 1D or 3D projection
"""
from __future__ import division
import math
import numpy

SQRT3_2 = math.sqrt(3) / 2


def draw_triangle_grid(ax, labels=['LAM', 'SAM', 'EAM'], linewidth=1.,
                       styles=['-', ':', '--'], fontsize=22):
    """
    .. versionchanged:: 0.2.0
        reworked internals and changed defaults

    Draws a triangle as axes, for a planar-simplex projection.

    Arguments:
        ax: an axis instance
        labels (iterable): list of string to be put for the axes
        styles (None, iterable): either *None* for solid lines or matplotlib
            line markers. These are in sync between the internal lines and
            the axes.
        linewidth (float): line width for the axes, the internal lines are
            equal to *0.75 * linewidth*
        fontsize (float): font size for the labels, the tick font size is
            equal to *0.75 * fontsize*
    """

    ax.set_ylim((-0.1, 1.0))
    ax.set_xlim((-0.1, 1.1))
    ax.set_axis_off()

    # SAM
    ax.text(0, -0.05, labels[0], fontsize=fontsize, ha='right')
    # LAM
    ax.text(1.0, -0.05, labels[1], fontsize=fontsize, ha='left')
    # EAB
    ax.text(0.5, SQRT3_2 + 0.025, labels[2], fontsize=fontsize, ha='center')

    values = numpy.arange(0.1, 1., 0.1)

    points = []

    # Bottom axis (LAM)
    x1, y1 = project_point((1, 0, 0))
    x2, y2 = project_point((0, 1, 0))
    points.append((
        [x1, x2],
        [y1, y2]
    ))
    # Left axis (EAM)
    x1, y1 = project_point((1, 0, 0))
    x2, y2 = project_point((0, 0, 1))
    points.append((
        [x1, x2],
        [y1, y2]
    ))
    # Right axis (SAM)
    x1, y1 = project_point((0, 1, 0))
    x2, y2 = project_point((0, 0, 1))
    points.append((
        [x1, x2],
        [y1, y2]
    ))

    if styles is None:
        styles = ['-'] * 3

    for (X, Y), style in zip(points, styles):
        ax.plot(X, Y, linestyle=style, linewidth=linewidth*1.5, color='k')

    for index in values:

        x1, y1 = project_point((index, 1 - index, 0))
        x2, y2 = project_point((0, 1 - index, index))
        # EAM
        ax.plot([x1, x2], [y1, y2], linestyle=styles[2], color='k',
                linewidth=linewidth)

        x3, y3 = project_point((0, index, 1 - index))
        x4, y4 = project_point((index, 0, 1 - index))
        # SAM
        ax.plot([x3, x4], [y3, y4], linestyle=styles[1], color='k',
                linewidth=linewidth)

        # LAM
        ax.plot([x4, x1], [y4, y1], linestyle=styles[0], color='k',
                linewidth=linewidth)

        # LAM labels
        ax.text(x1, y1 - 0.02, index, va='top', ha='center',
                fontsize=int(fontsize * 0.75))
        # EAM labels
        ax.text(x4 - 0.02, y4, 1 - index, ha='right', va='center',
                fontsize=int(fontsize * 0.75))
        # EAM labels
        ax.text(x3 + 0.02, y3, index, ha='left', va='center',
                fontsize=int(fontsize * 0.75))


def draw_1d_grid(ax, labels=['LAM', 'SAM'], fontsize=22):
    """
    .. versionchanged:: 0.2.0
        reworked internals and changed defaults

    Draws a 1D axis, to display propotions.

    Arguments:
        ax: an axis instance
        labels (iterable): list of string to be put for the axes
        fontsize (float): font size for the labels, the tick font size is
            equal to *0.75 * fontsize*
    """

    ax.set_ylim((-1.1, 1.1))
    ax.set_xlim((-0.1, 1.1))
    ax.set_axis_off()

    ax.text(-0.025, 0, labels[0], fontsize=fontsize, va='center', ha='right')
    ax.text(1.025, 0, labels[1], fontsize=fontsize, va='center', ha='left')

    ty = 1.0
    by = -1.0

    ax.plot([0, 1], [0, 0], linewidth=1, color='k')

    for x in numpy.arange(0.1, 1., 0.1):
        ax.plot([x, x], [by, ty], color='k', linestyle=':')
        ax.text(x, by - .1, x, fontsize=fontsize * 0.75, ha='center', va='top')
        ax.text(x, ty + .1, 1 - x, fontsize=fontsize * 0.75, ha='center',
                va='bottom')


def project_point(point):
    """
    Project a tuple containing coordinates (i.e. x, y, z) to planar-simplex.

    Arguments:
        point (tuple): contains the three coordinates to project

    Returns:
        tuple: the projected point in a planar-simplex
    """
    #a = point[0]
    b = point[1]
    c = point[2]

    x = b + (c / 2.)
    y = SQRT3_2 * c

    return (x, y)


def col_func_firstel(key, colors=None):
    if colors is None:
        return 'black'

    return colors[key[0]]


def col_func_name(key, func=None, colors=None):
    if (colors is None):
        return 'black'

    if func is not None:
        key = func(key)

    return colors[key]


def col_func_taxon(taxon_id, taxonomy, anc_ids, colpal):
    for anc_id, color in zip(anc_ids, colpal):
        if taxonomy.is_ancestor(taxon_id, anc_id):
            return color
    return 'black'


def draw_circles(ax, data, col_func=col_func_name, csize=200, alpha=.5,
                 sizescale=None, order=None, linewidths=0., edgecolor='none'):
    """
    .. versionchanged:: 0.2.0
        changed internals and added return value

    Draws a scatter plot over either a planar-simplex projection, if the number
    of coordinates is 3, or in a 1D axis.

    If the number of coordinates is 3, :func:`project_point` is used to project
    the point in 2 coordinates. The coordinates are converted in proportions
    internally.

    Arguments:
        ax: axis to plot on
        data (pandas.DataFrame): a DataFrame with 2 for a 1D plot or 3 columns
            for a planar-simplex
        col_func (func): a function that accept a parameter, an element of the
            DataFrame index and returns a colour for it
        csize (int): the base size of the circles
        alpha (float): transparency of the circles, between 0 and 1 included
        sizescale (None, pandas.Series): a Series or dictionary with the same
            elements as the Index of *data*, whose values are the size factors
            that are multiplied to *csize*. If **None**, the size of the
            circles is equal to *csize*
        order (None, iterable): iterable with the elements of *data* Index, to
            specify the order in which the circles must be plotted. If None,
            the order is the same as *data.index*
        linewidths (float): width of the circle line
        edgecolor (str): color of the circle line

    Returns:
        PathCollection: the return value of matplotlib *scatter*

    .. note::

        To **not** have circle lines, *edgecolor* must be *'none'* and
        *linewidths* equal *0*

    """
    #normalize data (ternary plots requires that the sum of a row is 1)
    data = data.div(data.sum(axis=1), axis=0)

    def no_project(a, b):
        if a > b:
            x = 1 - a
        else:
            x = b
        return x, 0.

    if order is None:
        order = data.index

    X = []
    Y = []
    for key, coord in data.loc[order].iterrows():
        # print coord
        x, y = project_point(coord) if len(coord) == 3 else no_project(*coord)
        X.append(x)
        Y.append(y)

    colors = [col_func(key) for key in order]

    paths = ax.scatter(
        X,
        Y,
        color=colors,
        alpha=alpha,
        linewidths=linewidths,
        edgecolor=edgecolor,
        marker=u'o',
        s=csize * (sizescale if sizescale is not None else 1)
    )
    return paths
