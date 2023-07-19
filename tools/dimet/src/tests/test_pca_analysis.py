from unittest import TestCase

import numpy as np

import pandas as pd

import processing.pca_analysis as pca_analysis


class TestPcaAnalysis(TestCase):

    def test_handle_nan_values_before_pca(self):
        data = {
            "c1": [0.1, 0, np.nan, 2],
            "c2": [25, 20, np.nan, 35],
            "c3": [np.nan, 5, np.nan, 4],
            "c4": [10, 12, np.nan, 34],
        }
        df = pd.DataFrame(data)
        result = pca_analysis.handle_nan_values_before_pca(df)
        self.assertEqual(result.shape[0], 3)
        self.assertEqual(result.shape[1], 4)
        self.assertTrue(
            np.isnan(df.loc[0, 'c3']) and
            result.loc[0, 'c3'] == 0.1
        )

    def test_reduce_data_df(self):
        data = {
            "c1": [0.1, 0.1, 2.0],
            "c2": [25.0, 20.0, 35.0],
            "c3": [0.1, 5.0, 4.0],
            "c4": [10.0, 12.0, 34.0]
        }
        df = pd.DataFrame(data)
        result = pca_analysis.reduce_data_df(df)
        self.assertTrue(
            np.allclose(np.array(result.loc[0, :]),
                        np.array([0.0098, 2.4536, 0.0098, 0.9814]), 4)
        )
        self.assertTrue(
            np.allclose(np.array(result.loc[1, :]),
                        np.array([0.0133, 2.6672, 0.6668, 1.6003]), 4)
        )
        self.assertTrue(
            np.allclose(np.array(result.loc[2, :]),
                        np.array([0.1268, 2.2194, 0.2536, 2.1560]), 4)
        )

    def test_compute_pca(self):
        data = {
            'beta-1': [0.0098, 0.0133, 0.1268],
            'beta-2': [2.4536, 2.6672, 2.2194],
            'ctrl-1': [0.0098, 0.6668, 0.2536],
            'ctrl-2': [0.9814, 1.6003, 2.1560]
        }
        df = pd.DataFrame(data)
        metadata_df = pd.DataFrame({
            'name_to_plot': ['beta-1', 'beta-2', 'ctrl-1', 'ctrl-2'],
            'condition': ['beta-glu', 'beta-glu', 'control', 'control']
        })
        pc_df, var_explained_df = pca_analysis.compute_pca(df, metadata_df)

        self.assertTrue(
            np.allclose(np.array(pc_df['PC1']),
                        np.array([-1.81, 2.34, -1.35, 0.82]), 2)
        )
        self.assertTrue(
            np.allclose(np.array(pc_df['PC2']),
                        np.array([-0.10, -0.38, -0.13, 0.62]), 2)
        )
        self.assertTrue(
            np.allclose(np.array(pc_df['PC3']),
                        np.array([-0.2279, -0.0277, 0.2562, -0.0004]), 4)
        )
        self.assertListEqual(list(pc_df['name_to_plot']),
                             list(metadata_df['name_to_plot']))
        self.assertListEqual(list(pc_df['condition']),
                             list(metadata_df['condition']))
        # var_explained_df :
        self.assertTrue(
            np.allclose(
                np.array(var_explained_df['Explained Variance %']),
                np.array([94.26, 4.47, 0.98]), 2
            )
        )
        self.assertListEqual(list(var_explained_df['PC']),
                             ["PC1", "PC2", "PC3"])
