#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

from unittest import TestCase

import helpers

import numpy as np

import pandas as pd


class TestHelpers(TestCase):
    def test_df_to_dict_bycomp(self):
        metadata_dict = {
            "name_to_plot": ["Ctrl_cell_0-1", "Ctrl_med_0-2",
                             "Ctrl_cell_0-2"],
            "condition": ["Control", "Control", "control"],
            "timepoint": ["T0", "T1", "T2"],
            "timenum": [0, 1, 2],
            "short_comp": ["cell", "medium", "cell"],
            "original_name": ["MCF001089_TD01", "MCF001089_TD02",
                              "MCF001089_TD03"],
        }
        abundances_dict = {
            "myindex": ["Fructose", "Fumaric_acid"],

            "MCF001089_TD01": [81.467, 1765.862],
            "MCF001089_TD02": [31.663, 2350.101],
            "MCF001089_TD03": [20.353, 19.404]
        }

        metadata_df = pd.DataFrame(metadata_dict)
        abundances_df = pd.DataFrame(abundances_dict)
        abundances_df = abundances_df.set_index(['myindex'])
        d = helpers.df_to_dict_by_compartment(df=abundances_df,
                                              metadata=metadata_df)
        self.assertEqual(list(d.keys()), ["cell", "medium"])
        self.assertEqual(d["cell"].shape, (2, 2))
        self.assertAlmostEqual(
            d["cell"].at['Fructose', "MCF001089_TD01"], 81.47, 2)

    def test_select_rows_by_fixed_values(self):
        data = {
            "Name": ["Alice", "Bob", "Charlie", "Dave", "Francoise"],
            "Age": [25, 30, 35, 40, 67],
            "City": ["London", "New York", "Paris", "London", "Paris"],
            "Country": ["UK", "USA", "France", "UK", "France"],
        }
        df = pd.DataFrame(data)

        # Specify the columns and values as lists
        columns_to_match = ["City", "Country"]
        values_to_match = [["London", "UK"], ["Paris", "France"]]

        # Select rows based on fixed values
        selected_rows = helpers.first_column_for_column_values(
            df, columns_to_match, values_to_match)

        self.assertEqual(selected_rows,
                         [["Alice", "Dave"], ["Charlie", "Francoise"]])

    def test_split_rows_by_threshold(self):
        data = {
            "Name": ["Alice", "Bob", "Charlie", "Dave", "Francoise"],
            "Age": [25, 30, 35, 40, 67],
            "City": ["London", "New York", "Paris", "London", "Paris"],
            "Country": ["UK", "USA", "France", "UK", "France"],
        }
        df = pd.DataFrame(data)
        df1, df2 = helpers.split_rows_by_threshold(df, "Age", 35)
        self.assertEqual(
            list(df1["Name"].values), ["Charlie", "Dave", "Francoise"]
        )  # assert might break due to ordering
        self.assertEqual(list(df2["Name"].values), ["Alice", "Bob"])

    def test_concatenate_dataframes(self):
        df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
        df2 = pd.DataFrame({"A": [10, 11, 12], "B": [13, 14, 15]})
        df3 = pd.DataFrame({"B": [16, 17, 18], "C": [19, 20, 21]})

        result = helpers.concatenate_dataframes(df1, df2, df3)
        result = result.fillna(-1)

        self.assertTrue(all(
            result["C"] == [7.0, 8.0, 9.0, -1.0, -1.0, -1.0, 19.0, 20.0,
                            21.0]))
        self.assertTrue(all(
            result["A"] == [1.0, 2.0, 3.0, 10.0, 11.0, 12.0, -1.0, -1.0,
                            0 - 1.0]))

    def test_row_wise_nanstd_reduction(self):
        data = {
            "A": [1, 0, 3, 4],
            "B": [5, 0, 6, 0],
            "C": [7, 0, 8, 9],
            "D": [10, 0, 11, 12],
        }
        df = pd.DataFrame(data)

        # Apply row-wise nanstd reduction
        result = helpers.row_wise_nanstd_reduction(df)
        self.assertTrue(np.allclose(np.array(result["B"]),
                                    np.array([1.529438, 0.0, 2.057983, 0.0])))
        self.assertTrue(np.allclose(np.array(result.iloc[1]),
                                    np.array([0.0, 0.0, 0.0, 0.0])))

    def test_compute_gmean_nonan(self):
        arr1 = np.array([1, 2, np.nan, 4, 5])
        arr2 = np.array([0, 0, 0, 0])
        gmean1 = helpers.compute_gmean_nonan(arr1)
        gmean2 = helpers.compute_gmean_nonan(arr2)
        self.assertAlmostEqual(gmean1, 2.514, 2)
        self.assertAlmostEqual(gmean2, np.finfo(float).eps)

    def test_first_column_for_column_values(self):
        data = {
            "Name": ["Alice", "Bob", "Charlie", "Dave"],
            "Age": [25, 30, 35, 40],
            "City": ["London", "Paris", "London", "New York"],
        }
        df = pd.DataFrame(data)

        # Define the columns and values to select
        columns = ["Age", "City"]
        values1 = [[30, "Paris"], [35, "London"]]
        values2 = [[30, "Paris", "France"], [35, "London"]]

        # Call the select_rows_by_fixed_values function
        result = helpers.first_column_for_column_values(df, columns, values1)
        self.assertEqual(result, [["Bob"], ["Charlie"]])
        self.assertRaises(ValueError, helpers.first_column_for_column_values,
                          df, columns, values2)

    def test_countnan_samples(self):
        data = {
            "c1": [0.1, np.nan, 2, np.nan],
            "c2": [25, np.nan, 35, 40],
            "c3": [np.nan, np.nan, np.nan, 365],
            "c4": [10, np.nan, 34, np.nan],
        }
        df = pd.DataFrame(data)
        result = helpers.countnan_samples(df, [["c1", "c2"], ["c3", "c4"]])
        self.assertTrue(np.any(
            np.array(result['count_nan_samples_group1']) == np.array(
                [0, 2, 0, 1])))
        self.assertTrue(np.any(
            np.array(result['count_nan_samples_group2']) == np.array(
                [1, 2, 1, 1])))

    def test_apply_multi_group_kruskal_wallis(self):
        np.random.seed(42)
        data = np.random.randint(0, 1000, size=(5, 12))
        df = pd.DataFrame(data,
                          columns=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
                                   'I', 'J', 'K', 'L'])

        cols = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I']]
        result = helpers.apply_multi_group_kruskal_wallis(df, cols)
        self.assertTrue(np.allclose(result['pvalue'].values,
                                    np.array([0.56, 0.56, 0.32, 0.07, 0.14]),
                                    2))

    def test_compute_padj(self):
        data = {
            'pvalue': [0.01, 0.02, np.nan, 0.03, 0.04],
            'other_column': [1, 2, 3, 4, 5]
        }
        df = pd.DataFrame(data)
        correction_alpha = 0.05
        correction_method = 'fdr_bh'
        df = helpers.compute_padj(df, correction_alpha,
                                  correction_method)
        self.assertAlmostEqual(df['padj'][1], 0.05, 3)
        self.assertTrue(np.isnan(df['padj'][2]))

    def test_verify_metadata_sample_not_duplicated(self):
        metadata = pd.DataFrame({
            "name_to_plot": ["Sample1", "Sample2", "Sample3", "Sample2",
                             "Sample4"]
        })
        self.assertRaises(ValueError,
                          helpers.verify_metadata_sample_not_duplicated,
                          metadata)
