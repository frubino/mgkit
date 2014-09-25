"""
Function to plot results from the datasets.

"""
from __future__ import division

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
import numpy
import logging

try:
    from scipy.spatial import ConvexHull
except ImportError:
    ConvexHull = None

try:
    import seaborn as sns
except ImportError:
    sns = None

LOG = logging.getLogger(__name__)

TAXON_COLOURS = {
    'archaea': '#4daf4a',
    'bacteria': '#377eb8',
    'fungi': '#ff7f00',
    'streptophyta': '#e41a1c',
    'ciliophora': '#984ea3'
}
"Default colours for root taxa"

#for backward compatibility
TAXON_COLORS = TAXON_COLOURS

DEFAULT_BOXPLOT_COLOURS = {
    'boxes': '#636363',
    'medians': '#f0f0f0',
    'whiskers': '#636363',
    'caps': 'black',
    'fliers': '#636363',
    'vals': '#636363',
}


DEFAULT_BOXPLOT_FONTCONF = {
    'rotation': 'vertical',
    'fontsize': 8
}


def lineplot_values_on_second_axis(gene_num, axis, colour='c', ylabel=''):
    """
    .. deprecated:: 0.1.13

    Adds a lineplot on a second axis using twinx
    """
    ax2 = axis.twinx()
    ax2.plot(axis.get_xticks(), gene_num, colour)

    ax2.set_ylabel(ylabel)
    ax2.set_yscale('log')

    return ax2


def plot_scatter_3d(data, labels, colours=None, pointsize=10., title='',
                    xlabel='', ylabel='', zlabel=''):
    """
    Scatter plot in 3d. Used for cluster results

    :param array data: :class:`numpy.array` with shape n, 3
    :param array labels: labels to categorise samples
    :param dict colours: dictionary whose keys are the labels and the values are
        valid :mod:`matplotlib` colours
    :param int pointsize: point size
    :param str title: plot title
    :param str xlabel: label for x axis
    :param str ylabel: label for y axis
    :param str zlabel: label for z axis

    :return: axis instance
    """

    fig = plt.figure()
    ax = plt.subplot(111, projection='3d')
    label_dict = {}
    for label, features in zip(labels, data):
        try:
            label_dict[label].append(features)
        except KeyError:
            label_dict[label] = [features]

    for label, data_points in label_dict.iteritems():
        if colours is None:
            colour = numpy.random.rand(3, 1)
        else:
            colour = colours[label]
        data_points = numpy.array(data_points)
        ax.scatter(data_points[:, 0], data_points[:, 1], data_points[:, 2],
                   s=pointsize, c=colour)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    fig.suptitle(title)

    return ax


def plot_scatter_2d(data, labels, colours=None, pointsize=10., title='',
                    xlabel='', ylabel='', centers=None, marker='*',
                    marker_colour=None, markersize=None, hull_points=True,
                    linewidth=0.2, legend=True, anno_center=True):
    """
    Scatter plot in 2d. Used for cluster results

    :param array data: :class:`numpy.array` with shape n, 2
    :param array labels: labels to categorise samples
    :param dict colours: dictionary whose keys are the labels and the values are
        valid :mod:`matplotlib` colours
    :param int pointsize: point size
    :param str title: plot title
    :param str xlabel: label for x axis
    :param str ylabel: label for y axis

    :return: axis instance
    """
    fig = plt.figure()
    ax = plt.subplot(111)
    label_dict = {}
    for label, features in zip(labels, data):
        try:
            label_dict[label].append(features)
        except KeyError:
            label_dict[label] = [features]

    if colours is None:
        #unify this, generate random colors in a go
        gen_colours = []

    colour_step = 256 // len(label_dict)

    for label, data_points in sorted(label_dict.iteritems()):
        if colours is None:
            colour = cm.Set1(label * colour_step)
            gen_colours.append(colour)
        else:
            colour = colours[label]
        data_points = numpy.array(data_points)
        if (len(data_points) > 2) and hull_points:
            if ConvexHull is None:
                LOG.error("ConvexHull from scipy.spatial couldn't be imported")
            else:
                hull = ConvexHull(data_points)
                for simplex in hull.simplices:
                    ax.plot(
                        data_points[simplex, 0],
                        data_points[simplex, 1],
                        '-',
                        colour=colour
                    )
        ax.scatter(
            data_points[:, 0],
            data_points[:, 1],
            s=pointsize,
            c=colour,
            linewidth=linewidth
        )

    if centers is not None:
        for idx, (center_x, center_y) in enumerate(centers):
            ax.plot(
                center_x,
                center_y,
                marker,
                markerfacecolor=marker_colour if marker_colour else gen_colours[idx],
                markeredgecolor='k',
                markersize=pointsize*1.2 if markersize is None else markersize,
                label=idx
            )
            if anno_center:
                ax.annotate(idx, xy=(center_x, center_y))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.suptitle(title)

    if legend:
        ax.legend()

    return ax


def barchart_categories(data, colours=None, title='', tickfont='small',
                        xlabel_dict=None, barlabel_dict=None, width=0.9,
                        rotation='vertical', file_name=None, fig_size=None,
                        fig_aspect=None):
    """
    :param data: DataFrame where the number of rows indicates how many bars will
        plotted per column
    :param colours: must be equal the number of data rows if supplied or it will
        be blue by default
    :param title: chart title
    :param tickfont: font size for ticks (only for column axis)
    :param xlabel_dict: a mapping to the acual labels to use for the columns.
        Defaults to columns' names
    :param barlabel_dict: a mapping to the acual labels to use for the row.
        Defaults to rows' names
    :param width: bar width

    :return: axis instance
    """

    if colours is None:
        colours = ('b',) * len(data.index)

    if fig_aspect is not None:
        fig_size = plt.figaspect(fig_aspect)
    fig = plt.figure(dpi=300, figsize=fig_size)

    ax = fig.add_subplot(111)

    base_index = numpy.arange(0.,
                              len(data.columns) * len(data.index),
                              len(data.index))

    for idx, (row_name, series) in enumerate(data.iterrows()):
        ax.bar(
            base_index + idx * width,
            series,
            color=colours[idx],
            label=row_name if barlabel_dict is None else barlabel_dict[idx]
        )

    ax.set_xticks(base_index + len(data.index) / 2.)

    xlabels = ax.set_xticklabels(
        [
            col_id if xlabel_dict is None else xlabel_dict[col_id]
            for col_id in data.columns
        ],
        rotation=rotation
    )

    ax.legend()

    ax.set_xlim(0, max(base_index) + (width * len(data.index)))
    ax.set_title(title)

    if not file_name is None:
        fig.savefig(file_name, bbox_inches='tight',
                    bbox_extra_artist=(xlabels,))
        fig.clf()

    return ax


def plot_contig_assignment_bar(series, taxon_colours=None, log_scale=False,
                               index=None, file_name=None, fig_aspect=None,
                               xlabels_size=8):
    """
    Plots barchart for contig assignment

    :param series: :class:`pandas.Series` instance with the data
    :param dict taxon_colours: colour of the bars for each taxon
    :param bool log_scale: if True the y axis is log scaled
    :param tuple fig_aspect: tuple with figure size
    :param int xlabels_size: size of the taxon labels
    :param index: optional :class:`pandas.Index` used to reindex the series
    :param str file_name: name of the file to write the graph to
    """
    LOG.info(
        "Plotting %d taxa, using colours (%s) and log scale (%s)",
        len(series),
        'no' if not taxon_colours else 'yes',
        'no' if not log_scale else 'yes'
    )

    if index is not None:
        series = series.reindex(index)

    if fig_aspect is None:
        fig_aspect = plt.figaspect(2. / len(series))

    LOG.info("Graph size %r and label size %d", tuple(fig_aspect), xlabels_size)

    fig = plt.figure(dpi=300, figsize=fig_aspect)

    ax = plt.subplot(111)

    colours = [taxon_colours[taxon_name] for taxon_name in series.index]

    ax.bar(range(1, series.count() + 1), series, color=colours,
           log=log_scale, align='center')

    ax.set_xlim(0, series.count() + 1)

    ax.set_xticks(range(1, series.count() + 1))
    xlabels = ax.set_xticklabels(series.index, rotation='vertical',
                                 fontsize=xlabels_size)
    ax.grid(True, axis='y')
    ax.set_ylabel('Number of contigs assigned', fontsize='medium',
                  horizontalalignment='right')
    ax.set_xlabel('Taxa', fontsize='medium')

    if file_name is not None:
        LOG.info("Saving graph to file %s", file_name)
        fig.savefig(file_name, bbox_inches='tight',
                    bbox_extra_artist=(xlabels,))
        fig.clf()

    return ax


def scatter_gene_values(gene_dict, xlabel="Profile pN/pS", ylabel="Rumen pN/pS",
                        title="", colours=None, file_name=None, plot_order=None,
                        line_colour='r', max_limit=None, axes=None):
    """
    Plots gene-taxon pN/pS from profiles against their observed values.

    :param dict gene_dict: dictionary that contains the data
    :param str xlabel: label for x axis
    :param str ylabel: label for y axis
    :param str title: graph title
    :param colours: colours used in for the different datasets; defaults to
        TAXON_COLOURS
    :param str file_name: path to which the graph is to be saved (by default)
        it doesn't write to disk
    :param iterable plot_order: the order used in plotting the data points;
        default to the order of the gene_dict dictionary keys
    :param coloUr: valid colour for the lines in the plot
    :param float max_limit: used to put a limit on the plot
    :param axes: optional axes used to draw the scatter plot

    :return: the axis object used for the plot
    """
    if colours is None:
        colours = TAXON_COLOURS

    if plot_order is None:
        plot_order = gene_dict.keys()

    if axes is None:
        fig = plt.figure()
        ax = plt.subplot(111)
    else:
        ax = axes

    for root in plot_order:
        root_values = gene_dict[root]
        values = numpy.array(root_values)

        if values.size == 0:
            continue

        ax.scatter(values[:, 0], values[:, 1], color=colours[root], label=root)

    if max_limit is None:
        max_limit = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot(numpy.arange(0., max_limit + 1), numpy.arange(0., max_limit + 1),
            color=line_colour)

    ax.axvline(1.0, color=line_colour)
    ax.axhline(1.0, color=line_colour)

    ax.set_xlim(left=0., right=max_limit)
    ax.set_ylim(bottom=0., top=max_limit)

    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if file_name is not None:
        LOG.info("Saving graph to file %s", file_name)
        fig.savefig(file_name)
        fig.clf()

    return ax


def boxplot_snp_dataframe(ratios, plot_order, taxon_colours=None,
                          labelfont='small', ylim=None, label_map=None,
                          file_name=None, fig_size=None,
                          title=None, log_scale=False, fig_aspect=None,
                          ylabel='', xlabel='', fig_multiplier=1.0,
                          fill_box=False, axes=None,
                          colours=None):
    """
    .. deprecated:: 0.1.7

    The function draws a series of boxplots from a DataFrame object, whose order
    is directed by the iterable plot_order. The columns of each DataFrame row
    contains the values for each boxplot. The function allows to write the plot
    to file and modify a few plot aspects. An axes object can be provided, which
    overrides fig_aspect, fig_size, fig_multiplier.

    :param ratios: dataframe to plot
    :param iterable plot_order: row order used to plot the boxes
    :param dict label_map: a map that converts the items in plot_order to a
        label used on the plot X axes
    :param str file_name: file name to save the figure to. Only used if axes
        is None
    :param str title: figure title
    :param str ylabel: label for the Y axes
    :param str xlabel: label fot the X axes
    :param str fill_box: if True the box is filled with the same colour of its
        outline

    :param int labelfont: font size for the x labels
    :param dict taxon_colours: a dictionary for each row index a colour is
        specified
    :param float ylim: limit for the Y axes
    :param fig_size: not used
    :param bool log_scale: if True the Y axes will be log scaled (not always
        possible)
    :param float fig_aspect: the aspect ratio of the figure axes, specified as a
        float (1.0 is a square)
    :param float fig_multiplier: it seems that the dpi is set when using
        figure.figaspect, so the actual inches used for the plot can be set as
        multiple of the results from figure.figaspect
    :param axes: a set of axes to be used. Nullifies fig_aspect, file_name,
        fig_size, fig_multiplier.
    :param str default_colour: default colour for boxes

    :return: a tuple with the axes used and the plot instance

    .. todo::

        Make the fig_aspect, fig_size, fig_multiplier values consistent.
    """

    if axes is None:
        if fig_aspect is not None:
            fig_size = plt.figaspect(fig_aspect) * fig_multiplier
        fig = plt.figure(figsize=fig_size)

        ax = plt.subplot(111)
    else:
        ax = axes

    plot = boxplot_dataframe(
        ratios,
        plot_order,
        ax,
        label_map=label_map,
        colours=colours,
        fill_box=fill_box,
        fonts={'rotation': 'vertical', 'fontsize': labelfont}
    )

    if not ylim is None:
        ax.set_ylim(0, ylim)

    if not title is None:
        ax.set_title(title)

    if log_scale:
        ax.set_yscale('log', basey=10)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if (file_name is not None) and (axes is None):
        LOG.info("Saving plot to file %s", file_name)
        fig.savefig(file_name, bbox_inches='tight',
                    bbox_extra_artist=(ax.get_xticklabels(),))
        fig.clf()

    return ax, plot


def boxplot_dataframe(dataframe, plot_order, axes, label_map=None, fonts=None,
                      fill_box=True, colours=None, data_colours=None,
                      box_vert=True):
    """
    .. versionadded:: 0.1.7
        To move from an all-in-one drawing to a more modular one.

    .. versionchanged:: 0.1.13
        added box_vert parameter

    The function draws a series of boxplots from a DataFrame object, whose order
    is directed by the iterable plot_order. The columns of each DataFrame row
    contains the values for each boxplot. An axes object is needed.

    :param dataframe: dataframe to plot
    :param iterable plot_order: row order used to plot the boxes
    :param axes: an axes instance
    :param dict label_map: a map that converts the items in plot_order to a
        label used on the plot X axes
    :param dict fonts: dictionary with properties for x axis labels,
        :data:`DEFAULT_BOXPLOT_FONTCONF` is used by default
    :param bool fill_box: if True each box is filled with the same colour of its
        outline
    :param dict colours: dictionary with properties for each boxplot if
        data_colours is None, whi overrides box, whiskers and fliers. Defaults
        to :data:`DEFAULT_BOXPLOT_COLOURS`
    :param dict data_colours: dictionary of colours for each boxplot, a set of
        colours can be obtained using func:`map_taxon_to_colours`
    :param bool box_vert: if False the boxplots are drawn horizontally

    :return: the plot data; same as matplotlib boxplot function
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

    plot_data = axes.boxplot(
        [dataframe.loc[x].dropna() for x in plot_order],
        vert=box_vert
    )

    for idx, row_id in enumerate(plot_order):
        box = plot_data['boxes'][idx]
        box.set_color(
            data_colours[row_id] if data_colours else colours['boxes']
        )
        if fill_box:
            box_coord = zip(box.get_xdata(), box.get_ydata())
            polygon = plt.Polygon(box_coord,
                                  facecolor=data_colours[row_id] if data_colours else colours['boxes'])
            axes.add_patch(polygon)

        plot_data['medians'][idx].set_color(colours['medians'])

    #It's got a different length (double the size of plot_order)
    for idx, tx in enumerate(plot_data['whiskers']):
        whisker = plot_data['whiskers'][idx]
        whisker.set_color(
            # data_colours[tx] if data_colours else colours['whiskers']
            colours['whiskers']
        )
        plot_data['caps'][idx].set_color(colours['caps'])
        flier = plot_data['fliers'][idx]
        flier.set_color(
            colours['whiskers']
            # data_colours[tx] if data_colours else colours['fliers']
        )

    if box_vert:
        ltick_setfunc = axes.set_xticklabels
        vtick_getfunc = axes.get_yticklabels
    else:
        ltick_setfunc = axes.set_yticklabels
        vtick_getfunc = axes.get_xticklabels

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


def map_taxon_to_colours(taxa, taxonomy, default_colour='#ffff33'):
    """
    Returns a dictionary of taxa and their assigned colours based on
    :data:`TAXON_COLOURS` and the taxonomy provided. Uses the
    :func:`taxon.UniprotTaxonomy.get_taxon_root` to determine the root of a
    taxon.

    :param iterable taxa: iterable of taxon ids
    :param UniprotTaxonomy taxonomy: taxonomy instance
    :param default_colour: colour used in case there's no known root for the
        taxon

    :return dict: dictionary mapping taxon_id to colour
    """
    LOG.info("Assigning taxa colours")
    tx_colours = {}

    for taxon_id in taxa:
        try:
            tx_colours[taxon_id] = TAXON_COLOURS[
                taxonomy.get_taxon_root(taxon_id).s_name
            ]
        except KeyError:
            tx_colours[taxon_id] = default_colour

    return tx_colours


get_taxon_colors_new = map_taxon_to_colours


def get_single_figure(dpi=300, figsize=(10, 20)):
    """
    Simple wrapper to init a single figure

    Arguments:
        dpi (int): dpi used for the figure
        figsize (tuple): size of the figure in inches

    Returns:
        tuple: the figure and axes objects
    """
    fig = plt.figure(dpi=dpi, figsize=figsize)
    ax = fig.add_subplot(111)
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


def boxplot_dataframe_multindex(dataframe, axes, plot_order=None, label_map=None,
                                fonts=None, fill_box=True, colours=None,
                                data_colours=None, box_vert=True):
    """
    .. versionadded:: 0.1.13

    .. todo::

        documentation

    The function draws a series of boxplots from a DataFrame object, whose order
    is directed by the iterable plot_order. The columns of each DataFrame row
    contains the values for each boxplot. An axes object is needed.

    :param dataframe: dataframe to plot
    :param iterable plot_order: row order used to plot the boxes
    :param axes: an axes instance
    :param dict label_map: a map that converts the items in plot_order to a
        label used on the plot X axes
    :param dict fonts: dictionary with properties for x axis labels,
        :data:`DEFAULT_BOXPLOT_FONTCONF` is used by default
    :param bool fill_box: if True each box is filled with the same colour of its
        outline
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
            box_coord = zip(box.get_xdata(), box.get_ydata())
            polygon = plt.Polygon(
                box_coord,
                facecolor=data_colours[category] if data_colours else colours['boxes']
            )
            axes.add_patch(polygon)

        plot_data['medians'][idx].set_color(colours['medians'])

    #It's got a different length (double the size of plot_order)
    for idx, tx in enumerate(plot_data['whiskers']):
        whisker = plot_data['whiskers'][idx]
        whisker.set_color(
            # data_colours[tx] if data_colours else colours['whiskers']
            colours['whiskers']
        )
        plot_data['caps'][idx].set_color(colours['caps'])
        flier = plot_data['fliers'][idx]
        flier.set_color(
            colours['whiskers']
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


def add_values_to_boxplot(dataframe, axes, plot_data, plot_order,
                          data_colours=None, alpha=0.5, s=80, marker='o',
                          linewidth=0.0):
    """
    .. versionadded:: 0.1.13

    Adds the values of a dataframe used in :func:`boxplot_dataframe` to the
    plot.

    Arguments:
        dataframe: dataframe with the values to plot
        axes: an axes instance
        plot_data: return value from :func:`boxplot_dataframe`
        plot_order (iterable): row order used to plot the boxes
        data_colours (dict): colors used for the values
        alpha (float): alpha value for the colour
        s (int): size of the marker drawn
        marker (str): one of the accepted matplotlib markers
        linewidth (float): width of the line used to draw the marker

    """
    for index, row_id in enumerate(plot_order):
        yvals = plot_data['medians'][index].get_ydata()
        mean_y = yvals.mean()

        axes.scatter(
            dataframe.loc[row_id].dropna(),
            [mean_y] * len(dataframe.loc[row_id].dropna()),
            c=DEFAULT_BOXPLOT_COLOURS['boxes'] if data_colours is None else data_colours[row_id],
            alpha=alpha,
            s=s,
            marker=marker,
            linewidth=linewidth,
            #this option put the dots below the lines of the boxplot
            zorder=1
        )
