import pandas as pd

import networkx as nx

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

            if tax_id not in taxid2data:
                taxid2data[tax_id] = {}

            if class_ == 'scientific name':
                taxid2data[tax_id]['scientific_name'] = name
            elif class_ == 'genbank common name':
                taxid2data[tax_id]['common_name'] = name

    # create network
    graph = nx.DiGraph()
    graph.name = 'Taxonomy'

    with open(fname_nodes) as fd:
        for line in tqdm(fd.readlines(), desc='Parsing nodes'):
            tax_id, parent_tax_id, rank, *_ = line.split('\t|\t')

            graph.add_node(tax_id, rank=rank, **taxid2data[tax_id])  # also updates if node exists
            graph.add_edge(parent_tax_id, tax_id)

    return graph


if __name__ == '__main__':
    graph = generate_taxonomy_network()

    for node in nx.shortest_path(graph.to_undirected(as_view=True), '9606', '4615'):
        print(node, graph.nodes[node])

    import IPython; IPython.embed()
