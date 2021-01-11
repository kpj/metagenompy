import pandas as pd
import pandas.testing as pdt

import pytest

import metagenompy


@pytest.mark.parametrize(
    'df_input,expected_output',
    [
        (
            pd.DataFrame(
                {
                    'taxid': ['1', '2', '2', '3', '3'],
                    'rank01': ['AA', 'AB', 'AB', 'BA', 'BA'],
                    'rank02': ['A', 'A', 'A', 'B', 'B'],
                }
            ),
            [
                pd.Series([3, 2], index=['A', 'B'], name='rank02'),
                pd.Series([2, 1, 2], index=['AB', 'AA', 'BA'], name='rank01'),
            ],
        ),
        (
            pd.DataFrame(
                {
                    'taxid': ['1', '1', '2'],
                    'rank01': ['AA', 'AA', pd.NA],
                    'rank02': ['A', 'A', 'B'],
                }
            ),
            [
                pd.Series([2, 1], index=['A', 'B'], name='rank02'),
                pd.Series([2, 1], index=['AA', pd.NA], name='rank01'),
            ],
        ),
    ],
)
def test_rank_frequency_computation(df_input, expected_output):
    freq_list = metagenompy.compute_rank_frequencies(
        df_input, rank_list=df_input.drop('taxid', axis=1).columns
    )

    assert len(freq_list) == len(expected_output)
    for freq, expected_freq in zip(freq_list, expected_output):
        pdt.assert_series_equal(freq, expected_freq)
