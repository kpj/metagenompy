import networkx as nx

import matplotlib.pyplot as plt


def plot_network(
    graph,
    ax=None,
    label_key='scientific_name',
    nodes_kws=dict(),
    labels_kws=dict(),
    edges_kws=dict(),
):
    """Visualize given taxonomy."""
    # compute additional properties
    node_labels = {
        n: data.get(label_key, '') for n, data in graph.nodes(data=True)
    }

    # create layout
    pos = nx.nx_agraph.pygraphviz_layout(graph, prog='dot')

    # visualize
    if ax is None:
        ax = plt.gca()

    nx.draw_networkx_nodes(graph, pos, ax=ax, **nodes_kws)
    nx.draw_networkx_labels(
        graph, pos, labels=node_labels, ax=ax, **labels_kws
    )
    nx.draw_networkx_edges(graph, pos, ax=ax, **edges_kws)

    ax.axis('off')
