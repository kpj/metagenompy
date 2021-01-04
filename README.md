# metagenompy

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
## 9606 {'rank': 'species', 'scientific_name': 'Homo sapiens', 'common_name': 'human'}
## 9605 {'rank': 'genus', 'scientific_name': 'Homo'}
## [..]
## 4614 {'rank': 'genus', 'scientific_name': 'Ananas'}
## 4615 {'rank': 'species', 'scientific_name': 'Ananas comosus', 'common_name': 'pineapple'}
```

### Easy transformation and visualization of taxonomy

`metagenompy`, e.g., allows you to quickly extract taxonomic entities of interest and visualize their relations.

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
metagenompy.plot_taxonomy(graph_zoom, ax=ax, labels_kws=dict(font_size=10))
fig.tight_layout()
fig.savefig('taxonomy.pdf')
```

![Exemplary taxonomy](gallery/taxonomy.png)
