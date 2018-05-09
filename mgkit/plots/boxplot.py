"""
.. versionadded:: 0.1.14

Code related to boxplots
"""
from __future__ import division
from builtins import zip
import logging
import numpy
from ..utils.common import deprecated
from .colors import float_to_hex_color

LOG = logging.getLogger(__name__)

from matplotlib.patches import Polygon

try:
    import seaborn as sns
except ImportError:
    sns = None

DEFAULT_BOXPLOT_FONTCONF = {
    'rotation': 'vertical',
    'fontsize': 8
}

DEFAULT_BOXPLOT_COLOURS = {
    'boxes': '#636363',
    'medians': '#f0f0f0',
    'whiskers': '#636363',
    'caps': 'black',
    'fliers': '#636363',
    'vals': '#636363',
}


def boxplot_dataframe_multindex(dataframe, axes, plot_order=None,
                                label_map=None, fonts=None, fill_box=True,
                                colours=None, data_colours=None,
                                box_vert=True):
    """
    .. versionadded:: 0.1.13

    .. todo::

        documentation

    The function draws a series of boxplots from a DataFrame object, whose
    order is directed by the iterable plot_order. The columns of each DataFrame
    row contains the values for each boxplot. An axes object is needed.

    :param dataframe: dataframe to plot
    :param iterable plot_order: row order used to plot the boxes
    :param axes: an axes instance
    :param dict label_map: a map that converts the items in plot_order to a
        label used on the plot X axes
    :param dict fonts: dictionary with properties for x axis labels,
        :data:`DEFAULT_BOXPLOT_FONTCONF` is used by default
    :param bool fill_box: if True each box is filled with the same colour of
        its outline
    :param dict colours: dictionary with properties for each boxplot if
        data_colours is None, whi overrides box, whiskers and fliers. Defaults
        to :data:`DEFAULT_BOXPLOT_COLOURS`
    :param dict data_colours: dictionary of colours for each boxplot, a set of
        colours can be obtained using func:`map_taxon_to_colours`

    :return: the plot data same as matplotlib boxplot function
    """

    if colours is not None:
        colours = dict(
            (feature, colours[feature]) if feature in colours else (feature, colour)
            for feature, colour in DEFAULT_BOXPLOT_COLOURS.items()
        )
        DEFAULT_BOXPLOT_COLOURS.copy().update(colours)
    else:
        colours = DEFAULT_BOXPLOT_COLOURS.copy()

    if fonts is not None:
        fonts = dict(
            (feature, fonts[feature]) if feature in fonts else (feature, option)
            for feature, option in DEFAULT_BOXPLOT_FONTCONF.items()
        )
        DEFAULT_BOXPLOT_FONTCONF.copy().update(fonts)
    else:
        fonts = DEFAULT_BOXPLOT_FONTCONF.copy()

    categories = set(dataframe.index.get_level_values(1))

    if (data_colours is None) and (sns is not None):
        data_colours = dict(
            zip(
                categories,
                sns.color_palette("hls", len(categories))
            )
        )

    if plot_order is None:
        plot_order = dataframe.index

    if label_map is None:
        label_map = []
        for label in dataframe.index.get_level_values(0):
            if label in label_map:
                continue
            label_map.append(label)

    plot_data = axes.boxplot(
        [dataframe.loc[x].dropna() for x in plot_order],
        vert=box_vert
    )

    for idx, row_id in enumerate(plot_order):
        category = row_id[1]
        box = plot_data['boxes'][idx]
        box.set_color(
            data_colours[category] if data_colours else colours['boxes']
        )
        if fill_box:
            box_coord = list(zip(box.get_xdata(), box.get_ydata()))
            polygon = Polygon(
                box_coord,
                facecolor=data_colours[category] if data_colours else colours['boxes']
            )
            axes.add_patch(polygon)

        plot_data['medians'][idx].set_color(colours['medians'])

    # It's got a different length (double the size of plot_order)
    for idx, tx in enumerate(plot_data['whiskers']):
        whisker = plot_data['whiskers'][idx]
        whisker.set_color(
            # data_colours[tx] if data_colours else colours['whiskers']
            colours['whiskers']
        )
        plot_data['caps'][idx].set_color(colours['caps'])

    for flier in plot_data['fliers']:
        flier.set_markerfacecolor(
            colours['fliers']
            # data_colours[tx] if data_colours else colours['fliers']
        )

    if box_vert:
        ltick_setfunc = axes.set_xticklabels
        vtick_getfunc = axes.get_yticklabels
        ptick_setfunc = axes.set_xticks
    else:
        ltick_setfunc = axes.set_yticklabels
        vtick_getfunc = axes.get_xticklabels
        ptick_setfunc = axes.set_yticks

    if fonts is not None:
        ptick_setfunc(
            numpy.arange(
                numpy.arange(1, len(categories) + 1).mean(),
                len(dataframe.index),
                len(categories)
            )
        )
        ltick_setfunc(
            label_map,
            rotation=fonts['rotation'],
            fontsize=fonts['fontsize']
        )
        for label in vtick_getfunc():
            label.set_fontsize(fonts['fontsize'])
    else:
        ltick_setfunc([])

    return plot_data


def add_values_to_boxplot(dataframe, ax, plot_data, plot_order,
                          data_colours=None, alpha=0.5, s=80, marker='o',
                          linewidth=0.01, box_vert=False):
    """
    .. versionadded:: 0.1.13

    .. versionchanged:: 0.1.14
        added *box_vert* parameter

    .. versionchanged:: 0.1.16
        changed default value for *linewidth*

    Adds the values of a dataframe used in :func:`boxplot_dataframe` to the
    plot. *linewidth* must be higher than 0 if a marker like *|* is used.

    A list of markers is available at
    `this page <http://matplotlib.org/api/markers_api.html>`_

    .. warning::

        Contrary to :func:`boxplot_dataframe`, the boxplot default is
        horizontal (*box_vert*). The default will change in a later version.

    Arguments:
        dataframe: dataframe with the values to plot
        ax: an axis instance
        plot_data: return value from :func:`boxplot_dataframe`
        plot_order (iterable): row order used to plot the boxes
        data_colours (dict): colors used for the values
        alpha (float): alpha value for the colour
        s (int): size of the marker drawn
        marker (str): one of the accepted matplotlib markers
        linewidth (float): width of the line used to draw the marker
        box_vert (bool): specify if the original boxplot is vertical or not
    """

    if data_colours is not None:
        # in case each color is a tuple of rgb floats, it converts them into
        # strings to avoid matplotlit to confuse the single color as different
        # shades of grey. It only happens when the number of data points in a
        # row is 3
        if not isinstance(data_colours[list(data_colours.keys())[0]], str):
            data_colours = dict(
                (key, float_to_hex_color(*value))
                for key, value in data_colours.items()
            )

    for index, row_id in enumerate(plot_order):
        if box_vert:
            xvals = plot_data['medians'][index].get_xdata()
            mean_x = xvals.mean()
            y = dataframe.loc[row_id].dropna()
            x = [mean_x] * dataframe.loc[row_id].count()
        else:
            yvals = plot_data['medians'][index].get_ydata()
            mean_y = yvals.mean()
            x = dataframe.loc[row_id].dropna()
            y = [mean_y] * dataframe.loc[row_id].count()

        if data_colours is None:
            DEFAULT_BOXPLOT_COLOURS['boxes']
        else:
            color = data_colours[row_id]

        ax.scatter(
            x,
            y,
            c=color,
            alpha=alpha,
            s=s,
            marker=marker,
            linewidth=linewidth,
            # this option put the dots below the lines of the boxplot
            zorder=1
        )


def boxplot_dataframe(dataframe, plot_order, ax, label_map=None, fonts=None,
                      fill_box=True, colours=None, data_colours=None,
                      box_vert=True, widths=0.5):
    """
    .. versionadded:: 0.1.7
        To move from an all-in-one drawing to a more modular one.

    .. versionchanged:: 0.1.13
        added box_vert parameter

    .. versionchanged:: 0.1.16
        added *widths* parameter

    The function draws a series of boxplots from a DataFrame object, whose
    order is directed by the iterable plot_order. The columns of each DataFrame
    row contains the values for each boxplot. An ax object is needed.

    :param dataframe: dataframe to plot
    :param iterable plot_order: row order used to plot the boxes
    :param ax: an axis instance
    :param dict label_map: a map that converts the items in plot_order to a
        label used on the plot X ax
    :param dict fonts: dictionary with properties for x axis labels,
        :data:`DEFAULT_BOXPLOT_FONTCONF` is used by default
    :param bool fill_box: if True each box is filled with the same colour of
        its outline
    :param dict colours: dictionary with properties for each boxplot if
        data_colours is None, whi overrides box, whiskers and fliers. Defaults
        to :data:`DEFAULT_BOXPLOT_COLOURS`
    :param dict data_colours: dictionary of colours for each boxplot, a set of
        colours can be obtained using func:`map_taxon_to_colours`
    :param bool box_vert: if False the boxplots are drawn horizontally
    :param float widths: width (scalar or array) of the boxplots width(s)

    :return: the plot data; same as matplotlib boxplot function
    """

    if colours is not None:
        colours = dict(
            (feature, colours[feature]) if feature in colours else (feature, colour)
            for feature, colour in DEFAULT_BOXPLOT_COLOURS.items()
        )
        # DEFAULT_BOXPLOT_COLOURS.copy().update(colours)
    else:
        colours = DEFAULT_BOXPLOT_COLOURS.copy()

    if fonts is not None:
        fonts = dict(
            (feature, fonts[feature]) if feature in fonts else (feature, option)
            for feature, option in DEFAULT_BOXPLOT_FONTCONF.items()
        )
        DEFAULT_BOXPLOT_FONTCONF.copy().update(fonts)
    else:
        fonts = DEFAULT_BOXPLOT_FONTCONF.copy()

    plot_data = ax.boxplot(
        [dataframe.loc[x].dropna() for x in plot_order],
        vert=box_vert,
        widths=widths
    )

    for idx, row_id in enumerate(plot_order):
        box = plot_data['boxes'][idx]
        box.set_color(
            data_colours[row_id] if data_colours else colours['boxes']
        )
        if fill_box:
            box_coord = list(zip(box.get_xdata(), box.get_ydata()))
            polygon = Polygon(
                box_coord,
                facecolor=data_colours[row_id] if data_colours else colours['boxes']
            )
            ax.add_patch(polygon)

        plot_data['medians'][idx].set_color(colours['medians'])

    # It's got a different length (double the size of plot_order)
    for idx, tx in enumerate(plot_data['whiskers']):
        whisker = plot_data['whiskers'][idx]
        whisker.set_color(
            # data_colours[tx] if data_colours else colours['whiskers']
            colours['whiskers']
        )
        plot_data['caps'][idx].set_color(colours['caps'])

    for flier in plot_data['fliers']:
        flier.set_color(
            colours['whiskers']
            # data_colours[tx] if data_colours else colours['fliers']
        )

    if box_vert:
        ltick_setfunc = ax.set_xticklabels
        vtick_getfunc = ax.get_yticklabels
    else:
        ltick_setfunc = ax.set_yticklabels
        vtick_getfunc = ax.get_xticklabels

    if fonts is not None:
        ltick_setfunc(
            [
                label if label_map is None else label_map[label]
                for label in plot_order
            ],
            rotation=fonts['rotation'], fontsize=fonts['fontsize']
        )
        for label in vtick_getfunc():
            label.set_fontsize(fonts['fontsize'])
    else:
        ltick_setfunc([])

    return plot_data


def add_significance_to_boxplot(sign_indices, ax, pos, box_vert=True,
                                fontsize=16):
    """
    .. versionadded:: 0.1.16

    Add significance groups to boxplots

    Arguments:
        sign_indices (iterable): iterable in which each element is a tuple;
            each element of the tuple is the numerical index of the position of
            the significant boxplot
        ax: an axis instance
        pos (tuple): the 2 values are the coordinates for the top line, and the
            the lowest bound for the whisker
        box_vert (bool): if the boxplot is vertical
        fontsize (float): size for the * (star)
    """
    maxval, spine = pos
    text = maxval

    for index, (sign1, sign2) in enumerate(sign_indices):
        spacer = (index * ((maxval - spine) * 2.5))
        factors = [spine, maxval]
        x = numpy.array(
            factors + factors[::-1]
        ) - spacer

        y = [sign1 + 1, sign1 + 1, sign2 + 1, sign2 + 1]

        if box_vert:
            x, y = y, x
        ax.plot(x, y, linestyle='-', c='k', alpha=.75)
        if box_vert:
            xtext, ytext = numpy.mean(x), text - spacer
        else:
            xtext, ytext = text - spacer, numpy.mean(y)
        ax.text(xtext, ytext, '*', fontsize=fontsize)

__all__ = [
    'add_values_to_boxplot',
    'add_significance_to_boxplot',
    'boxplot_dataframe_multindex',
    'boxplot_dataframe'
]
