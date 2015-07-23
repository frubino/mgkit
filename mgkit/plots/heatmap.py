"""
.. versionadded:: 0.1.14

Code related to heatmaps.
"""
import logging
import numpy
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as hclust
from .utils import get_grid_figure

LOG = logging.getLogger(__name__)


def baseheatmap(data, ax, norm=None, cmap=None, xticks=None, yticks=None,
                fontsize=18, meshopts=None):
    """
    A basic heatmap using :func:`matplotlib.pyplot.pcolormesh`. It expect a
    :class:`pandas.DataFrame`.

    .. note::

        Rows a plot bottom to up, while the columns left to right. Change the
        order of the DataFrame if needed.

    Arguments:
        data (pandas.DataFrame): matrix to plot. The DataFrame labels are used
        ax: axes to use
        norm: if needed, :class:`matplotlib.colors.BoundaryNorm` or
            :class:`matplotlib.colors.Normalize` can be used to fine tune the
            colors
        cmap (None, matplotlib.colors.ListedColormap): color map to use
        xticks (None, dict): dictionary with additional options to pass to
            *set_xticklabels*
        yticks (None, dict): dictionary with additional options to pass to
            *set_yticklabels*
        fontsize (int): font size to use for the labels
        meshopts (None, dict): additional options to pass to
            :func:`matplotlib.pyplot.pcolormesh`

    Returns:
        matplotlib.collections.QuadMesh: the return value of
        :func:`matplotlib.pyplot.pcolormesh`
    """
    mesh_args = {
        'edgecolor': 'face'
    }
    xticks_args = {
        'pos': 'bottom',
        'rotation': 'vertical'
    }
    yticks_args = {
        'pos': 'right',
        'rotation': 'horizontal'
    }

    if meshopts is not None:
        mesh_args.update(meshopts)

    if xticks is not None:
        xticks_args.update(xticks)

    if yticks is not None:
        yticks_args.update(yticks)

    mesh = ax.pcolormesh(data.values, cmap=cmap, norm=norm, **mesh_args)

    ax.xaxis.set_ticks_position(xticks_args['pos'])
    ax.set_xticks(numpy.arange(0.5, len(data.columns), 1))
    ax.set_xticklabels(
        data.columns, fontsize=fontsize, rotation=xticks_args['rotation']
    )

    ax.yaxis.set_ticks_position(yticks_args['pos'])
    ax.set_yticks(numpy.arange(0.5, len(data.index), 1))
    ax.set_yticklabels(
        data.index, fontsize=fontsize, rotation=yticks_args['rotation']
    )

    ax.set_ylim(top=len(data))
    ax.set_xlim(right=len(data.columns))

    return mesh


def grouped_spine(groups, labels, ax, which='y', spine='right',
                  spine_opts=None):
    """
    .. versionchanged:: 0.2.0
        added *va*, *ha* keys to *spine_opts*, changed the label positioning

    Changes the spine of an heatmap axis given the groups of labels.

    .. note::

        It should work for any plot, but was not tested

    Arguments:
        groups (iterable): a nested list where each is element is a list
            containing the labels belong to that group.
        labels (iterable): an iterable with the labels of the groups. Needs to
            be in the same order as groups
        ax: axis to use (same as heatmap)
        which (str): to specify the axis, either *x* or *y*
        spine (str): position of the spine. if *which* is **x** accepted values
            are *top* and *bottom*, if which is **y** *left* and *right* are
            accepted
        spine_opts (dict): additional options to pass to the spine class

    """
    spine_args = dict(
        ec='k',
        position=15,
        group_lw=1.5,
        in_length=13,
        out_length=15,
        fontsize=18,
        va='center',
        ha='left'
    )

    if spine_opts is not None:
        spine_args.update(spine_opts)

    major_ticks = [0]
    minor_ticks = []

    for group in groups:
        if not group:
            continue

        major_ticks.append(major_ticks[-1] + len(group))

        group_half = len(group) // 2
        addendum = 0.5

        if len(group) % 2 == 0:
            addendum = 0.
        minor_ticks.append(
            group_half + major_ticks[-2] + addendum
        )

    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_position(('outward', spine_args['position']))
    ax.spines[spine].set_edgecolor(spine_args['ec'])
    ax.spines[spine].set_linewidth(spine_args['group_lw'])

    if which is 'y':
        axis = ax.yaxis
    else:
        axis = ax.xaxis

    axis.set_ticks(major_ticks)
    axis.set_ticks(minor_ticks, minor=True)
    axis.set_ticks_position(spine)
    axis.set_ticklabels(labels, minor=True, va=spine_args['va'],
                        ha=spine_args['ha'])
    axis.set_ticklabels([], minor=False)
    axis.set_tick_params(
        direction='in',
        length=spine_args['in_length'],
        which='major',
        width=spine_args['group_lw']
    )
    axis.set_tick_params(
        direction='out',
        length=spine_args['out_length'],
        which='minor',
        width=spine_args['group_lw'],
        pad=spine_args['fontsize'],
        labelsize=spine_args['fontsize']
    )


def dendrogram(data, ax, method='complete', orientation='top', use_dist=True,
               dist_func=distance.pdist):
    """
    .. versionchanged:: 0.1.16
        added *use_dist* and *dist_func* parameters

    Plots a dendrogram of the clustered rows of the given matrix; if the
    columns are to be clustered, the transposed matrix needs to be passed.

    Arguments:
        data (pandas.DataFrame): matrix to plot. The DataFrame labels are used
        ax: axes to use
        method (str): clustering method used, internally
            :func:`scipy.cluster.hierarchy.linkage` is used.
        orientation (str): direction for the plot. *top*, *bottom*, *left* and
            *right* are accepted; *top* will draw the leaves at the bottom.
        use_dist (bool): if True, the function *dist_func* will be applied to
            *data* to get a distance matrix
        dist_func (func): distance function to be used

    Returns:
        The dendrogram plotted, as returned by
        :func:`scipy.cluster.hierarchy.dendrogram`

    """
    if use_dist:
        data = dist_func(data)

    pairwise_dists = distance.squareform(data)

    clusters = hclust.linkage(pairwise_dists, method=method)

    dendrogram = hclust.dendrogram(
        clusters,
        ax=ax,
        link_color_func=lambda x: 'black',
        orientation=orientation
    )
    ax.grid(False)
    ax.set_axis_bgcolor('white')
    ax.set_xticks([])
    ax.set_yticks([])

    return dendrogram


def heatmap_clustered(data, figsize=(10, 5), cmap=None, norm=None):
    """
    Plots a heatmap clustered on both rows and columns.

    Arguments:
        data (pandas.DataFrame): matrix to plot. The DataFrame labels are used
        figsize (tuple): passed to :func:`mgkit.plots.utils.get_grid_figure`
        cmap (None, matplotlib.colors.ListedColormap): color map to use
        norm: if needed, :class:`matplotlib.colors.BoundaryNorm` or
            :class:`matplotlib.colors.Normalize` can be used to fine tune the
            colors
    """

    fig, gs = get_grid_figure(
        2,
        2,
        figsize=figsize,
        wspace=0,
        hspace=0,
        width_ratios=[0.25, 1],
        height_ratios=[0.25, 1]
    )
    dendr = dendrogram(data, fig.add_subplot(gs[1, 0]), orientation='right')
    dendc = dendrogram(data.T, fig.add_subplot(gs[0, 1]))

    baseheatmap(
        data.ix[dendr['leaves'], dendc['leaves']],
        fig.add_subplot(gs[1, 1]),
        cmap=cmap,
        norm=norm,
        xticks={'pos': 'bottom'}
    )


__all__ = ['baseheatmap', 'grouped_spine', 'dendrogram', 'heatmap_clustered']
