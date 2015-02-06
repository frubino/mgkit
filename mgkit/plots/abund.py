"""
.. versionadded:: 0.1.15

Module to plot relative abundances in a 1D or 3D projection
"""
from __future__ import division
from shapely import geometry
import math
import numpy

SQRT3_2 = math.sqrt(3) / 2


def draw_triangle_grid(ax, labels=['LAB', 'SAB', 'EAB'],
                       styles=['-.', ':', '--'], linewidth=1., fontsize=16):
    """
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

    lines = geometry.LineString(
        [(0, 0), (0.5, SQRT3_2), (1, 0), (0, 0)]
    )

    if styles is None:
        # solid lines for the axes
        ax.plot(*lines.xy, linewidth=linewidth, color='k')
    else:
        for index, style in zip(range(3), styles):
            ax.plot(
                lines.xy[0][index:index+2],
                lines.xy[1][index:index+2],
                linewidth=linewidth,
                color='k',
                linestyle=style
            )

    ax.set_ylim((-0.1, 1.0))
    ax.set_xlim((-0.1, 1.1))
    ax.set_axis_off()

    # font size for the labels
    fontdict = {'fontsize': fontsize}

    ax.text(-0.05, -0.05, labels[0], fontdict=fontdict)
    ax.text(1.0, -0.05, labels[1], fontdict=fontdict)
    ax.text(0.5 - 0.025, SQRT3_2 + 0.05, labels[2], fontdict=fontdict)

    # font size for the ticks
    fontdict = {'fontsize': fontsize * 0.75}

    Xvalues = numpy.arange(0.1, SQRT3_2, SQRT3_2 / 10)
    Yvalues = numpy.arange(0.1, 1., 0.1)

    for x, y in zip(Xvalues, Yvalues):
        line = geometry.LineString([(0, x), (1, x)])
        points = line.intersection(lines)
        line = geometry.LineString(points)
        ax.plot(
            *line.xy,
            color='k',
            linestyle=':' if styles is None else styles[0],
            linewidth=linewidth * 0.75
        )
        ax.text(
            points.geoms[0].x - 0.05,
            points.geoms[0].y,
            "{0:.0f}".format(y * 100),
            fontdict=fontdict
        )

        line = geometry.LineString([points.geoms[0], (y, 0)])
        ax.plot(
            *line.xy,
            color='k',
            linestyle=':' if styles is None else styles[2],
            linewidth=linewidth * 0.75
        )
        ax.text(y, - 0.05, "{0:.0f}".format((1 - y) * 100), fontdict=fontdict)

        line = geometry.LineString([points.geoms[1], (1.-y, 0)])
        ax.plot(
            *line.xy,
            color='k',
            linestyle=':' if styles is None else styles[1],
            linewidth=linewidth * 0.75
        )
        ax.text(
            points.geoms[1].x + 0.025,
            points.geoms[1].y,
            "{0:.0f}".format((1 - y) * 100),
            fontdict=fontdict
        )


def draw_1d_grid(ax, labels=['LAB', 'SAB'], fontsize=16):
    """
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

    # fontsize for labels
    fontdict = {'fontsize': fontsize}

    ax.text(-0.075, -.05, labels[0], fontdict=fontdict)
    ax.text(1.025, -.05, labels[1], fontdict=fontdict)

    # fontsize for ticks
    fontdict = {'fontsize': fontsize * 0.75}

    ty = 1.0
    by = -1.0

    lines = geometry.LineString([(0, 0), (1, 0)])
    ax.plot(*lines.xy, linewidth=1, color='k')

    for x in numpy.arange(0.1, 1., 0.1):
        line = geometry.LineString([(x, by), (x, ty)])
        ax.plot(*line.xy, color='k', linestyle=':')
        ax.text(
            x - 0.01,
            by - .75,
            "{0:.0f}".format(x * 100),
            fontdict=fontdict
        )
        ax.text(
            x - 0.01,
            ty + .25,
            "{0:.0f}".format((1 - x) * 100),
            fontdict=fontdict
        )


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

    .. note::

        To **not** have circle lines, *edgecolor* must be *'none'* and
        *linewidths* equal *0*

    """
    #normalize data (ternary plots requires that the sum of a row is 1)
    data = data.div(data.sum(axis=1), axis=0)

    if order is None:
        order = data.index

    for key in order:
        coord = data.loc[key]
        c = col_func(key)
        #ternary plot
        if len(coord) == 3:
            x, y = project_point(coord)
        #1D plot
        else:
            a, b = coord
            if a > b:
                x = 1 - a
            else:
                x = b
            y = 0.
        ax.scatter(
            x,
            y,
            c=c,
            alpha=alpha,
            linewidths=linewidths,
            edgecolor=edgecolor,
            marker=u'o',
            s=csize * (sizescale[key] if sizescale is not None else 1)
        )
