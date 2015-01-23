"""
.. versionadded:: 0.1.12

Graph module
"""

import itertools
import networkx as nx


def build_graph(id_links, name, edge_type='', weight=0.5):
    """
    .. versionadded:: 0.1.12

    Builds a networkx graph from a dictionary of nodes, as outputted by
    :meth:`mgkit.kegg.KeggClientRest.get_pathway_links`. The graph is
    undirected, and all edges weight are the same.

    Arguments:
        id_links (dict): dictionary with the links
        name (str): name of the graph
        edge_type (str): an optional name for the `edge_type` attribute
            set for each edge
        weight (float): the weight assigned to each edge in the graph

    Returns:
        graph: an instance of :class:`networkx.Graph`
    """
    g = nx.Graph()
    g.name = name

    for id1, id2s in id_links.iteritems():
        g.add_node(id1, id=id1)
        for id2 in id2s:
            g.add_node(id2, id=id2)
            g.add_edge(id1, id2, edge_type=edge_type, weight=weight)

    return g


def build_weighted_graph(id_links, name, weights, edge_type=''):
    """
    .. versionadded:: 0.1.14

    Builds a networkx graph from a dictionary of nodes, as outputted by
    :meth:`mgkit.kegg.KeggClientRest.get_pathway_links`. The graph is
    undirected, and all edges weight are the same.

    Arguments:
        id_links (dict): dictionary with the links
        name (str): name of the graph
        edge_type (str): an optional name for the `edge_type` attribute
            set for each edge
        weight (float): the weight assigned to each edge in the graph

    Returns:
        graph: an instance of :class:`networkx.Graph`
    """
    g = nx.Graph()
    g.name = name

    for id1, id2s in id_links:
        g.add_node(id1, id=id1)
        for id2 in id2s:
            g.add_node(id2, id=id2)
            try:
                weight = weights[(id1, id2)]
            except KeyError:
                weight = weights.get(id1, 0.5)
            g.add_edge(id1, id2, edge_type=edge_type, weight=float(weight))

    return g


def copy_nodes(g, graph1):
    """
    .. versionadded:: 0.1.12

    Used by :func:`link_nodes` to copy nodes
    """
    for node, data in graph1.nodes_iter(data=True):
        g.add_node(
            "{0}_{1}".format(graph1.name, data['id']),
            module=graph1.name,
            alpha=1.0,
            **data
        )


def copy_edges(g, graph1):
    """
    .. versionadded:: 0.1.12

    Used by :func:`link_nodes` to copy edges
    """
    for node1, node2, data in graph1.edges_iter(data=True):
        g.add_edge(
            "{0}_{1}".format(graph1.name, node1),
            "{0}_{1}".format(graph1.name, node2),
            data
        )


def link_nodes(g, graph1, graph2, id_filter, link_type, weight):
    """
    .. versionadded:: 0.1.12

    Used by :func:`link_graph` to link nodes with the same *id*
    """
    nodes1 = set(
        data['id']
        for node, data in graph1.nodes_iter(data=True)
        if id_filter(data['id'])
    )
    nodes2 = set(
        data['id']
        for node, data in graph2.nodes_iter(data=True)
        if id_filter(data['id'])
    )

    for node_id in (nodes1 & nodes2):
        g.node["{0}_{1}".format(graph1.name, node_id)]['alpha'] = 0.5
        g.node["{0}_{1}".format(graph2.name, node_id)]['alpha'] = 0.5
        g.add_edge(
            "{0}_{1}".format(graph1.name, node_id),
            "{0}_{1}".format(graph2.name, node_id),
            edge_type=link_type,
            weight=weight
        )

EDGE_LINKS = [
    (lambda x: x.startswith('C'), 'CP_LINK', 0.0),
    (lambda x: x.startswith('K'), 'KO_LINK', 0.0)
]
"Sample edge_links for :func:`link_graph`"


def link_graph(graphs, edge_links):
    """
    .. versionadded:: 0.1.12

    Link nodes of a set of graphs using the specifics in edge_links.
    The resulting graph nodes are renamed, and the nodes that are shared
    between the graphs linked.

    Arguments:
        graphs: iterable of graphs
        edge_links: iterable with function, edge_type and weight for the
            links between graphs

    Returns:
        graph: an instance of :class:`networkx.Graph`
    """
    g = nx.Graph()

    for graph in graphs:
        copy_nodes(g, graph)
        copy_edges(g, graph)

    for graph1, graph2 in itertools.combinations(graphs, r=2):
        for id_filter, link_type, weight in edge_links:
            link_nodes(g, graph1, graph2, id_filter, link_type, weight)

    return g


def filter_graph(graph, id_list, filter_func=lambda x: x.startswith('K')):
    """
    .. versionadded:: 0.1.12

    Filter a graph based on the `id_list` provided and the filter function
    used to test the id attribute of each node.

    A node is removed if `filter_func` returns True on a node and its id
    attribute is not in `id_list`

    Arguments:
        graph: the graph to filter
        id_list (iterable): the list of nodes that are to remain in the
            graph
        filter_func (func): function which accept a single parameter and
            return a boolean

    Returns:
        graph: an instance of :class:`networkx.Graph`
    """
    graph = graph.copy()
    nodes = [
        node
        for node, data in graph.nodes_iter(data=True)
        if filter_func(data['id']) and (data['id'] not in id_list)
    ]

    graph.remove_nodes_from(nodes)

    graph.remove_nodes_from(nx.isolates(graph))

    return graph


def annotate_graph_nodes(graph, attr, id_map, default=None):
    """
    .. versionadded:: 0.1.12

    Add/Changes nodes attribute `attr` using a dictionary of ids->values.

    .. note::

        If the id is not found in `id_map`:

        * default is None: no value added for that node
        * default is not None: the node attribute will be set to `default`

    Arguments:
        graph: the graph to annotate
        attr (str): the attribute to annotate
        id_map (dict): the dictionary with the values for each node
        default: the value used in case an `id` is not found in `id_map`
    """
    for node, data in graph.nodes_iter(data=True):
        try:
            data[attr] = id_map[data['id']]
        except KeyError:
            if default is not None:
                data[attr] = default
