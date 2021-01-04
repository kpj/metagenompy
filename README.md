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
>>> import metagenompy
>>> import networkx as nx

>>> graph = metagenompy.generate_taxonomy_network()
>>> print(nx.info(graph))
Name: Taxonomy
Type: DiGraph
Number of nodes: 2298250
Number of edges: 2298250
Average in degree:   1.0000
Average out degree:   1.0000
```
