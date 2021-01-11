import importlib.resources as pkg_resources

import pandas as pd

from . import example_data


def load_example_dataset(name='blast_result.csv.gz'):
    with pkg_resources.path(example_data, name) as path:
        return pd.read_csv(path)
