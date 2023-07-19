#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
import os
from typing import List, Tuple, Union

from constants import (
    assert_literal,
    data_files_keys_type,
)

from data import Dataset

import helpers

import numpy as np

from omegaconf import DictConfig

import pandas as pd

from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)


def handle_nan_values_before_pca(df: pd.DataFrame) -> pd.DataFrame:
    """
    1. drops rows having NaN in all values
    2. imputes NaN remaining values with the min of the dataframe
    """
    df = df.dropna(how="all", axis=0)
    df = df.fillna(df[df > 0].min().min())
    return df


def reduce_data_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Row-wise reduction of quantitative dataframe for pca
    """
    df_red = helpers.row_wise_nanstd_reduction(df)  # reduce rows
    if np.isinf(np.array(df_red)).any():  # avoid Inf error
        return df
    else:
        return df_red


def compute_pca(
        quantitative_df: pd.DataFrame,
        metadata_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Uses sklearn PCA.
    Returns the two computed metrics:
     - Computed dimensions coeficients (pc_df)
     - Explained Variance in percentages (var_explained_df)
    """
    dims = min(metadata_df.shape[0],
               quantitative_df.shape[0])  # min(nb_samples, nb_features)

    X = np.transpose(np.array(quantitative_df))
    pca = PCA(n_components=dims)
    pc = pca.fit_transform(X)
    pc_df = pd.DataFrame(data=pc,
                         columns=['PC' + str(i) for i in range(1, dims + 1)])
    pc_df = pc_df.assign(name_to_plot=quantitative_df.columns)
    pc_df = pd.merge(pc_df, metadata_df, on='name_to_plot')

    var_explained_df = pd.DataFrame({
        'Explained Variance %': pca.explained_variance_ratio_ * 100,
        'PC': ['PC' + str(i) for i in range(1, dims + 1)]})

    return pc_df, var_explained_df


def pca_on_split_dataset(compartment_df: pd.DataFrame,
                         metadata_co_df: pd.DataFrame, chosen_column: str,
                         description: List[str]):
    """
    Using compartment specific dataframes,
    splits metadata and data by selected column ("condition" or "timepoint"),
    and computes PCA on each subset.
    The results are added to the dictionary of results.
    """
    assert len(metadata_co_df['short_comp'].unique()) == 1
    assert chosen_column in ["condition", "timepoint"]
    unique_nominal_values = metadata_co_df[chosen_column].unique().tolist()
    pca_tables_dict = {}
    for cond_or_timepoint in unique_nominal_values:
        # example of cond_or_timepoint : 'T0'
        metadata_co_ct_df = metadata_co_df.loc[
            metadata_co_df[chosen_column] == cond_or_timepoint, :]
        df = compartment_df[metadata_co_ct_df['name_to_plot']]
        df = handle_nan_values_before_pca(df)
        df = reduce_data_df(df)
        pc_df, var_explained_df = compute_pca(df, metadata_co_ct_df)

        key_unique = tuple(
            [description[0],  # index 0 is datatype (eg: abundances)
             cond_or_timepoint,
             description[1]]  # index 1 is compartment
        )
        pca_tables_dict[key_unique] = {
              'pc': pc_df,
              'var': var_explained_df
        }

    return pca_tables_dict


def pca_global_compartment_dataset(df: pd.DataFrame,
                                   metadata_co_df: pd.DataFrame,
                                   description: List[str]):
    df = df[metadata_co_df['name_to_plot']]
    df = handle_nan_values_before_pca(df)
    df = reduce_data_df(df)
    pc_df, var_explained_df = compute_pca(df, metadata_co_df)

    pca_tables_dict = {tuple([description[0], description[1]]): {
        'pc': pc_df,
        'var': var_explained_df
    }}
    return pca_tables_dict


def send_to_tables(pca_results_compartment_dict: dict,
                   out_table_dir: str) -> None:
    """ Save each result to .csv files """
    for tup in pca_results_compartment_dict.keys():
        out_table = "--".join(list(tup))
        for df in pca_results_compartment_dict[tup].keys():
            pca_results_compartment_dict[tup][df].to_csv(
                os.path.join(out_table_dir, f"{out_table}_{df}.csv"),
                sep='\t', index=False)
    logger.info(f"Saved pca tables in {out_table_dir}")


def run_pca_analysis(file_name: data_files_keys_type,
                     dataset: Dataset, cfg: DictConfig,
                     out_table_dir: str, mode: str) -> Union[None, dict]:
    """
    Generates all PCA results, both global (default) and with splited data.
     - mode='save_tables', the PCA tables are saved to .csv;
     or
     - mode='return_results_dict', returns the results object (dict)
    """
    assert_literal(file_name, data_files_keys_type, "file name")

    metadata_df = dataset.metadata_df

    pca_results_dict = dict()

    impute_value = cfg.analysis.method.impute_values[file_name]
    for compartment, compartmentalized_df in \
            dataset.compartmentalized_dfs[file_name].items():
        df = compartmentalized_df
        df = df[(df.T != 0).any()]  # delete rows being zero all values
        val_instead_zero = helpers.arg_repl_zero2value(impute_value, df)
        df = df.replace(to_replace=0, value=val_instead_zero)

        metadata_co_df = metadata_df[metadata_df['short_comp'] == compartment]

        pca_compartment_dict = pca_global_compartment_dataset(
            df, metadata_co_df, description=[file_name, compartment]
        )
        pca_results_dict.update(pca_compartment_dict)

        if cfg.analysis.method.pca_split_further is not None:
            for column in cfg.analysis.method.pca_split_further:
                pca_results_split_data_dict = pca_on_split_dataset(
                     df, metadata_co_df,
                     column,  description=[file_name, compartment]
                )
                pca_results_dict.update(pca_results_split_data_dict)
            # end for
        # end if
    # end for

    if mode == "save_tables":
        send_to_tables(pca_results_dict, out_table_dir)

    elif mode == "return_results_dict":
        return pca_results_dict
