"""
.. versionadded:: 0.1.12

Graph module
"""

import itertools
from future.utils import viewitems
from xml.etree import ElementTree
from . import DependencyError

try:
    import networkx as nx
except ImportError:
    raise DependencyError('networkx')


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

    for id1, id2s in viewitems(id_links):
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


def rename_graph_nodes(graph, name_func=None, exclude_ids=None):
    new_graph = nx.Graph()

    if exclude_ids is None:
        exclude_ids = set()

    for node, data in graph.nodes_iter(data=True):
        new_graph.add_node(
            node if node in exclude_ids else name_func(node),
            **data
        )

    for node1, node2, data in graph.edges_iter(data=True):
        new_graph.add_edge(
            node1 if node1 in exclude_ids else name_func(node1),
            node2 if node2 in exclude_ids else name_func(node2),
            **data
        )

    return new_graph


def copy_nodes(g, graph1, name=None, id_attr=None, **kwd):
    """
    .. versionadded:: 0.1.12

    Used by :func:`link_nodes` to copy nodes
    """

    if name is None:
        name = graph1.name

    if id_attr is None:
        id_attr = 'id'

    for node, data in graph1.nodes_iter(data=True):
        data = data.copy()
        data.update(kwd)
        g.add_node(
            "{0}_{1}".format(name, data[id_attr]),
            **data
        )


def copy_edges(g, graph1, name=None, **kwd):
    """
    .. versionadded:: 0.1.12

    Used by :func:`link_nodes` to copy edges
    """

    if name is None:
        name = graph1.name

    for node1, node2, data in graph1.edges_iter(data=True):
        data = data.copy()
        data.update(kwd)
        g.add_edge(
            "{0}_{1}".format(name, node1),
            "{0}_{1}".format(name, node2),
            **data
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


def from_kgml(entry, graph=None, rn_ids=None):
    """
    .. versionadded:: 0.3.1

    Given a KGML file (as string), representing a pathway in Kegg, returns a
    networkx DiGraph, using reaction directionality included in the KGML. If a
    reaction is reversible, 2 edges (from and to) for each compound/reaction
    pair are added, giving the bidirectionality.

    .. note::

        substrate and products included in a KGML don't represent the complete
        reaction, excluding in general cofactors or more general terms.
        Those can be added using :func:`add_module_compounds`, which may be
        more useful when used with a restricted number of reactions (e.g.
        a module)

    Arguments:
        entry (str): KGML file as a string, or anything that can be passed to
            ElementTree
        graph (graph): an instance of a networkx DiGraph if the network is to
            be updated with a new KGML, if `None` a new one is created
        rn_ids (set): a set/list of reaction IDs that are to be included, if
            `None` all reactions are used

    Returns:
        graph: a networkx DiGraph with the reaction/compounds
    """
    if graph is None:
        graph = nx.DiGraph()

    for entry in ElementTree.fromstring(entry).findall('reaction'):
        # the "reaction" defined is equivalent to a EC number, meaning multiple
        # reactions IDs in kegg may belong to it. They are stored in the name
        # attribute, separated by space
        reactions = entry.attrib['name'].split(' ')

        for reaction in reactions:
            # Each reaction ID is (as usual) prepended by the type "rn", which
            # we don't need
            reaction = reaction.split(':')[1]
            if (rn_ids is not None) and (reaction not in rn_ids):
                continue
            # definition of the reaction direction, by manual either reversible
            # or irreversible
            if entry.attrib['type'] == 'irreversible':
                reversible = False
            else:
                reversible = True

            substrates = []
            products = []

            graph.add_node(
                reaction,
                node_type='reaction',
                reaction_type='reversible' if reversible else 'irreversible'
            )

            # separating substrate and products from the compounds listed in
            # the reaction definition
            for compound in entry:
                cpd_id = compound.attrib['name'].split(':')[1]
                if compound.tag == 'substrate':
                    substrates.append(cpd_id)
                else:
                    products.append(cpd_id)

            for substrate in substrates:
                if substrate not in graph:
                    graph.add_node(substrate, node_type='substrate')
                else:
                    # cases where the substrate is product of other reactions
                    if graph.node[substrate]['node_type'] != 'substrate':
                        graph.node[substrate]['node_type'] = 'mixed'

                graph.add_edge(
                    substrate,
                    reaction,
                    reaction_type='reversible' if reversible else 'irreversible'
                )
                # if reversible add the reciprocal edge
                if reversible:
                    graph.add_edge(
                        reaction,
                        substrate,
                        reaction_type='reversible' if reversible else 'irreversible'
                    )

            for product in products:
                if product not in graph:
                    graph.add_node(product, node_type='product')
                else:
                    # cases where the product is product of other reactions
                    if graph.node[product]['node_type'] != 'product':
                        graph.node[product]['node_type'] = 'mixed'

                graph.add_edge(
                    reaction,
                    product,
                    reaction_type='reversible' if reversible else 'irreversible'
                )
                # if reversible add the reciprocal edge
                if reversible:
                    graph.add_edge(
                        product,
                        reaction,
                        reaction_type='reversible' if reversible else 'irreversible'
                    )

    return graph


def add_module_compounds(graph, rn_defs):
    """
    .. versionadded:: 0.3.1

    Modify in-place a graph, by adding additional compounds from a dictionary
    of definitions. It uses the reversible/irreversible information for each
    reaction to add the correct number of edges to the graph.

    Arguments:
        graph (graph): a graph to update with additional compounds
        rn_defs (dict): a dictionary, whose keys are reactions IDs and the
            values are instances of :class:`mgkit.kegg.KeggReaction`
    """
    for rn_id, (left_cp, right_cp) in viewitems(rn_defs):

        reaction_type = graph.node[rn_id]['reaction_type']
        reversible = True if reaction_type == 'reversible' else False
        for cp_id in left_cp | right_cp:
            if ((rn_id, cp_id) in graph.edges()) or ((rn_id, cp_id) in graph.edges()):
                continue
            if cp_id not in graph:
                graph.add_node(cp_id, node_type='addition')
            if reversible:
                graph.add_edge(rn_id, cp_id, reaction_type=reaction_type)
                graph.add_edge(cp_id, rn_id, reaction_type=reaction_type)
            else:
                if cp_id in graph.predecessors(rn_id):
                    graph.add_edge(cp_id, rn_id, reaction_type=reaction_type)
                else:
                    graph.add_edge(rn_id, cp_id, reaction_type=reaction_type)
