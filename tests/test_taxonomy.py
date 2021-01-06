import pandas as pd

import pytest

import metagenompy


@pytest.fixture
def taxdump(tmp_path):
    fname_nodes = tmp_path / 'nodes.dmp'
    with open(fname_nodes, 'w') as fd:
        # tax_id, parent_tax_id, rank, more stuff...
        fd.write("""1\t|\t1\t|\tno rank\t|\tfoo\t|\n""")
        fd.write("""2\t|\t1\t|\tkingdom\t|\tfoo\t|\n""")
        fd.write("""3\t|\t1\t|\tkingdom\t|\tfoo\t|\n""")
        fd.write("""20\t|\t2\t|\tphylum\t|\tfoo\t|\n""")
        fd.write("""30\t|\t3\t|\tclade\t|\tfoo\t|\n""")
        fd.write("""300\t|\t30\t|\tphylum\t|\tfoo\t|\n""")

    fname_names = tmp_path / 'names.dmp'
    with open(fname_names, 'w') as fd:
        # tax_id, name, unique_name, name_class
        fd.write("""1\t|\troot\t|\t\t|\tcommon name\t|\n""")
        fd.write("""2\t|\tfoo2\t|\t\t|\tcommon name\t|\n""")
        fd.write("""3\t|\tfoo3\t|\t\t|\tcommon name\t|\n""")
        fd.write("""20\t|\tbar20\t|\t\t|\tcommon name\t|\n""")
        fd.write("""30\t|\tbar30\t|\t\t|\tcommon name\t|\n""")
        fd.write("""300\t|\tbaz300\t|\t\t|\tcommon name\t|\n""")
        fd.write("""300\t|\thui\t|\t\t|\tspecial\t|\n""")
        fd.write("""300\t|\tbuh\t|\t\t|\tspecial\t|\n""")

    return fname_nodes, fname_names


def test_taxonomy_creation(taxdump):
    fname_nodes, fname_names = taxdump
    graph = metagenompy.generate_taxonomy_network(
        fname_nodes=fname_nodes, fname_names=fname_names
    )

    expected_nodes = {
        '1': {'rank': 'no rank', 'common_name': 'root'},
        '2': {'rank': 'kingdom', 'common_name': 'foo2'},
        '3': {'rank': 'kingdom', 'common_name': 'foo3'},
        '20': {'rank': 'phylum', 'common_name': 'bar20'},
        '30': {'rank': 'clade', 'common_name': 'bar30'},
        '300': {
            'rank': 'phylum',
            'common_name': 'baz300',
            'special': ['hui', 'buh'],
        },
    }
    expected_edges = [
        ('1', '1', {}),
        ('1', '2', {}),
        ('1', '3', {}),
        ('2', '20', {}),
        ('3', '30', {}),
        ('30', '300', {}),
    ]

    assert dict(graph.nodes(data=True)) == expected_nodes
    assert list(graph.edges(data=True)) == expected_edges


def test_taxonomy_condensation(taxdump):
    fname_nodes, fname_names = taxdump
    graph = metagenompy.generate_taxonomy_network(
        fname_nodes=fname_nodes, fname_names=fname_names
    )

    metagenompy.condense_taxonomy(graph)

    expected_nodes = {
        '1': {'rank': 'no rank', 'common_name': 'root'},
        '2': {'rank': 'kingdom', 'common_name': 'foo2'},
        '3': {'rank': 'kingdom', 'common_name': 'foo3'},
        '20': {'rank': 'phylum', 'common_name': 'bar20'},
        '300': {
            'rank': 'phylum',
            'common_name': 'baz300',
            'special': ['hui', 'buh'],
        },
    }
    expected_edges = [
        ('1', '1', {}),
        ('1', '2', {}),
        ('1', '3', {}),
        ('2', '20', {}),
        ('3', '300', {}),
    ]

    assert dict(graph.nodes(data=True)) == expected_nodes
    assert list(graph.edges(data=True)) == expected_edges


def test_node_highlighting(taxdump):
    fname_nodes, fname_names = taxdump
    graph = metagenompy.generate_taxonomy_network(
        fname_nodes=fname_nodes, fname_names=fname_names
    )

    graph_sub = metagenompy.highlight_nodes(graph, ['2', '30'])

    expected_nodes = {
        '1': {'rank': 'no rank', 'common_name': 'root'},
        '2': {'rank': 'kingdom', 'common_name': 'foo2'},
        '3': {'rank': 'kingdom', 'common_name': 'foo3'},
        '30': {'rank': 'clade', 'common_name': 'bar30'},
    }
    expected_edges = [
        ('1', '1', {}),
        ('1', '2', {}),
        ('1', '3', {}),
        ('3', '30', {}),
    ]

    assert dict(graph_sub.nodes(data=True)) == expected_nodes
    assert list(graph_sub.edges(data=True)) == expected_edges


def test_classification(taxdump):
    fname_nodes, fname_names = taxdump
    graph = metagenompy.generate_taxonomy_network(
        fname_nodes=fname_nodes, fname_names=fname_names
    )

    assert metagenompy.classify_taxid(graph, '2', 'kingdom') == '2'

    assert metagenompy.classify_taxid(graph, '20', 'kingdom') == '2'
    assert metagenompy.classify_taxid(graph, '30', 'kingdom') == '3'

    assert pd.isna(metagenompy.classify_taxid(graph, '30', 'phylum'))
