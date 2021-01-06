import itertools

import pandas as pd
import networkx as nx

from tqdm import tqdm


tqdm.pandas()


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

            graph.add_node(
                tax_id, rank=rank, **taxid2data[tax_id]
            )  # also updates if node exists
            graph.add_edge(parent_tax_id, tax_id)

    return graph


def condense_taxonomy(
    graph,
    abbreviated_lineage=[
        'no rank',
        'superkingdom',
        'kingdom',
        'phylum',
        'class',
        'order',
        'suborder',
        'family',
        'genus',
        'species',
        'subspecies',
    ],
):
    """Contract certain edges in-place."""
    node_data = list(graph.nodes(data=True))
    for node, data in tqdm(node_data):
        if data['rank'] not in abbreviated_lineage:
            sources = graph.predecessors(node)
            targets = graph.successors(node)

            new_edges = [
                (source, target)
                for source, target in itertools.product(sources, targets)
                if source != target
            ]

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
            return pd.NA

        node = parents[0]


def classify_dataframe(
    graph,
    df,
    rank_list=['species', 'genus', 'class', 'superkingdom'],
    name_key='scientific_name',
):
    """Create dataframe with various rank classifications."""

    def func(taxid, rank):
        clf_id = classify_taxid(graph, taxid, rank)

        if pd.isna(clf_id) or name_key is None:
            return clf_id

        return graph.nodes[clf_id][name_key]

    for rank in tqdm(rank_list, desc='Classifying'):
        df[rank] = df['taxid'].apply(func, rank=rank)

    return df


def aggregate_classifications(
    df, rank, min_fraction=1, group_column='qseqid', taxon_column='taxid'
):
    """Aggregate dataframe at given rank."""
    assert min_fraction > 0.5, 'min_fraction must be greater than 0.5'

    def func(x):
        # discard entry if there is no majority
        freqs = x[rank].value_counts(normalize=True)
        if freqs.max() < min_fraction:
            # returning pd.NA leads to "AttributeError: 'NAType' object has no attribute 'index'"
            return None

        # retain only most common taxon
        top_taxon = x[taxon_column].value_counts().index[0]
        return x[x[taxon_column] == top_taxon].iloc[0]

    return (
        df.groupby(group_column)
        .progress_apply(func)
        .fillna(pd.NA)
        .drop(group_column, axis=1)
    )


def highlight_nodes(graph, node_list, root_node='1'):
    """Retain only specified nodes and their paths to root."""
    node_subset = set()
    for node in node_list:
        path = nx.shortest_path(graph, root_node, node)
        node_subset.update(path)

    return graph.subgraph(node_subset)
