# metagenompy

[![PyPI](https://img.shields.io/pypi/v/metagenompy.svg?style=flat)](https://pypi.python.org/pypi/metagenompy)
[![Tests](https://github.com/kpj/metagenompy/workflows/Tests/badge.svg)](https://github.com/kpj/metagenompy/actions)

Your all-inclusive package for aggregating and visualizing metagenomic BLAST results.


## Installation

```bash
$ pip install metagenompy
```


## Usage

### NCBI taxonomy as NetworkX object

The core of `metagenompy` is a taxonomy as a networkX object.
This means that all your favorite algorithms work right out of the box.

```python
import metagenompy
import networkx as nx


# load taxonomy
graph = metagenompy.generate_taxonomy_network()

# print path from human to pineapple
for node in nx.shortest_path(graph.to_undirected(as_view=True), '9606', '4615'):
    print(node, graph.nodes[node])
## 9606 {'rank': 'species', 'authority': 'Homo sapiens Linnaeus, 1758', 'scientific_name': 'Homo sapiens', 'genbank_common_name': 'human', 'common_name': 'man'}
## 9605 {'rank': 'genus', 'authority': 'Homo Linnaeus, 1758', 'scientific_name': 'Homo', 'common_name': 'humans'}
## [..]
## 4614 {'rank': 'genus', 'authority': 'Ananas Mill., 1754', 'scientific_name': 'Ananas'}
## 4615 {'rank': 'species', 'authority': ['Ananas comosus (L.) Merr., 1917', 'Ananas lucidus Mill., 1754'], 'scientific_name': 'Ananas comosus', 'synonym': ['Ananas comosus var. comosus', 'Ananas lucidus'], 'genbank_common_name': 'pineapple'}
```

### Easy transformation and visualization of taxonomy

Extract taxonomic entities of interest and visualize their relations:

```python
import metagenompy
import matplotlib.pyplot as plt


# load and condense taxonomy to relevant ranks
graph = metagenompy.generate_taxonomy_network()
metagenompy.condense_taxonomy(graph)

# highlight interesting nodes
graph_zoom = metagenompy.highlight_nodes(graph, [
    '9606',  # human
    '9685',  # cat
    '9615',  # dog
    '4615',  # pineapple
    '3747',  # strawberry
    '4113',  # potato
])

# visualize result
fig, ax = plt.subplots(figsize=(10, 10))
metagenompy.plot_network(graph_zoom, ax=ax, labels_kws=dict(font_size=10))
fig.tight_layout()
fig.savefig('taxonomy.pdf')
```

<img src="gallery/taxonomy.png" width="50%">


Classify taxonomic entities at different ranks:

```python
import metagenompy
import pandas as pd


# load taxonomy
graph = metagenompy.generate_taxonomy_network()

# classification
tmp = []
for taxid in ['9606', '9685', '3747']:
    for rank in ['class', 'order']:
        clf_id = metagenompy.classify_taxid(graph, taxid, rank)
        tmp.append({
            'taxid': graph.nodes[taxid]['scientific_name'],
            'rank': rank,
            'clf': graph.nodes[clf_id]['scientific_name']
        })

pd.DataFrame(tmp)
##                  taxid   rank            clf
## 0         Homo sapiens  class       Mammalia
## 1         Homo sapiens  order       Primates
## 2          Felis catus  class       Mammalia
## 3          Felis catus  order      Carnivora
## 4  Fragaria x ananassa  class  Magnoliopsida
## 5  Fragaria x ananassa  order        Rosales
```
