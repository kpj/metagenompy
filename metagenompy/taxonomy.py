import itertools

import pandas as pd

import networkx as nx

import matplotlib.pyplot as plt

from tqdm import tqdm


def generate_taxonomy_network(
    fname_nodes='nodes.dmp', fname_names='names.dmp'
):
    """Create taxonomic network."""
    # associate tax ids with additional data
    taxid2data = {}
    with open(fname_names) as fd:
        for line in tqdm(fd.readlines(), desc='Parsing names'):
            tax_id, name, _, class_ = line.split('\t|\t')
            class_ = class_[:-3]  # remove '\t|\n' suffix
            class_ = class_.replace(' ', '_')

            if tax_id not in taxid2data:
                taxid2data[tax_id] = {}

            # create list of multiple names exist for same class of same taxid
            if class_ in taxid2data[tax_id]:
                if isinstance(taxid2data[tax_id][class_], list):
                    taxid2data[tax_id][class_].append(name)
                else:
                    tmp = taxid2data[tax_id][class_]
                    taxid2data[tax_id][class_] = [tmp, name]
            else:
                taxid2data[tax_id][class_] = name

    # create network
    graph = nx.DiGraph()
    graph.name = 'Taxonomy'

    with open(fname_nodes) as fd:
        for line in tqdm(fd.readlines(), desc='Parsing nodes'):
            tax_id, parent_tax_id, rank, *_ = line.split('\t|\t')

            graph.add_node(tax_id, rank=rank, **taxid2data[tax_id])  # also updates if node exists
            graph.add_edge(parent_tax_id, tax_id)

    return graph


def condense_taxonomy(
    graph,
    relevant_ranks=['no rank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']
):
    """Contract certain edges in-place."""
    node_data = list(graph.nodes(data=True))
    for node, data in tqdm(node_data):
        if data['rank'] not in relevant_ranks:
            sources = graph.predecessors(node)
            targets = graph.successors(node)

            new_edges = [(source, target)
                         for source, target in itertools.product(sources, targets)
                         if source != target]

            graph.add_edges_from(new_edges)
            graph.remove_node(node)


def classify_taxid(graph, taxid, rank):
    """Determine name of given taxid at specified rank."""
    node = taxid
    while True:
        current_rank = graph.nodes[node]['rank']
        if current_rank == rank:
            return node

        parents = list(graph.predecessors(node))
        assert len(parents) == 1

        if node == parents[0]:
            raise RuntimeError(f'Cannot classify {taxid} at rank {rank}')

        node = parents[0]


def highlight_nodes(graph, node_list, root_node='1'):
    """Retain only specified nodes and their paths to root."""
    node_subset = set()
    for node in node_list:
        path = nx.shortest_path(graph, root_node, node)
        node_subset.update(path)

    return graph.subgraph(node_subset)


def plot_taxonomy(
    graph,
    ax=None, label_key='scientific_name',
    nodes_kws=dict(), labels_kws=dict(), edges_kws=dict()
):
    """Visualize given taxonomy."""
    # compute additional properties
    node_labels = {n: data.get(label_key, '')
                   for n, data in graph.nodes(data=True)}

    # create layout
    pos = nx.nx_agraph.pygraphviz_layout(graph, prog='dot')

    # visualize
    if ax is None:
        ax = plt.gca()

    nx.draw_networkx_nodes(graph, pos, ax=ax, **nodes_kws)
    nx.draw_networkx_labels(graph, pos, labels=node_labels, ax=ax, **labels_kws)
    nx.draw_networkx_edges(graph, pos, ax=ax, **edges_kws)

    ax.axis('off')
