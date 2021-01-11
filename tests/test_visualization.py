import pandas as pd
import pandas.testing as pdt

import metagenompy


def test_rank_frequency_computation():
    df = pd.DataFrame(
        {
            'taxid': ['1', '2', '2', '3', '3'],
            'rank01': ['AA', 'AB', 'AB', 'BA', 'BA'],
            'rank02': ['A', 'A', 'A', 'B', 'B'],
        }
    )
    freq_dict = metagenompy.compute_rank_frequencies(
        df, rank_list=['rank01', 'rank02']
    )

    assert list(freq_dict.keys()) == ['rank02', 'rank01']
    pdt.assert_series_equal(
        freq_dict['rank01'],
        pd.Series([2, 1, 2], index=['AB', 'AA', 'BA'], name='rank01'),
    )
    pdt.assert_series_equal(
        freq_dict['rank02'], pd.Series([3, 2], index=['A', 'B'], name='rank02')
    )
