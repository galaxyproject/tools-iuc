#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
import logging
from collections.abc import Iterable
from functools import reduce
from typing import Dict, List

from constants import assert_literal, overlap_methods_types

import numpy as np

import pandas as pd

import scipy
from scipy import stats

import statsmodels.stats.multitest as ssm

logger = logging.getLogger(__name__)


def compute_padj(df: pd.DataFrame, correction_alpha: float,
                 correction_method: str) -> pd.DataFrame:
    '''
    Performs multiple hypothesis testing correction on the p-values in the
    DataFrame.
    Deals with the situation where pvalue column can contain np.nan values
    Adds a new column called "padj" with the adjusted p-values.
    '''
    tmp = df.copy()
    # inspired from R documentation in p.adjust :
    tmp["pvalue"] = tmp[["pvalue"]].fillna(1)

    (sgs, corrP, _, _) = ssm.multipletests(tmp["pvalue"],
                                           alpha=float(correction_alpha),
                                           method=correction_method)
    df = df.assign(padj=corrP)
    truepadj = []
    for v, w in zip(df["pvalue"], df["padj"]):
        if np.isnan(v):
            truepadj.append(v)
        else:
            truepadj.append(w)
    df = df.assign(padj=truepadj)

    return df


def row_wise_nanstd_reduction(df: pd.DataFrame) -> pd.DataFrame:
    """
    Performs a row-wise reduction of the DataFrame by dividing each value
    by its row's standard deviation,
    considering the presence of NaN values.
    """
    std_values = df.apply(lambda row: np.nanstd(row), axis=1)
    std_values[
        std_values == 0] = 1
    # Replace zero standard deviations with 1 to avoid division by zero
    result = df.div(std_values, axis=0)
    return result


def concatenate_dataframes(df1: pd.DataFrame, df2: pd.DataFrame,
                           df3: pd.DataFrame) -> pd.DataFrame:
    """
    Concatenate df2 and df2 to df1 ; fill the missing values with np.nan
    """
    assert set(df2.columns).issubset(set(df1.columns))
    assert set(df3.columns).issubset(set(df1.columns))
    df2 = df2.reindex(columns=df1.columns, fill_value=np.nan)
    df3 = df3.reindex(columns=df1.columns, fill_value=np.nan)
    # please leave ignore_index as False:
    # otherwise numbers and not metabolites appear in .csv exported results:
    result = pd.concat([df1, df2, df3], ignore_index=False)
    return result


def split_rows_by_threshold(df: pd.DataFrame, column_name: str,
                            threshold: float) -> (pd.DataFrame, pd.DataFrame):
    """
    Splits the dataframe into rows having column_name value >= threshold and
     the rest
    Returns two dataframes
    """
    message = "Error in split_rows_by_threshold \
    check qualityDistanceOverSpan parameter in the analysis YAML file"

    try:
        good_df = df.loc[df[column_name] >= threshold, :]
        undesired_rows = set(df.index) - set(good_df.index)
        bad_df = df.loc[list(undesired_rows)]
    except Exception as e:
        logger.info(e)
        logger.info(message)

    return good_df, bad_df


def calculate_gmean(df: pd.DataFrame,
                    groups: List[List[str]]) -> pd.DataFrame:
    """
    Calculates the geometric mean for each row in the specified
     column groups and adds the corresponding values
    as new columns to the DataFrame. Additionally, adds a column
    with the ratio of the two geometric means.
    Takes care of the potential division by zero error
    by replacing 0 by 1e-10 in division

    groups: A list with two sublists containing column names from df.

    Returns:
        The modified DataFrame with additional columns for the calculated
        geometric means
        and the ratio of means (FC - Fold Change).
    """
    for i, group in enumerate(groups):
        gmean_col = f"gmean_{i + 1}"
        df[gmean_col] = df[group].apply(lambda x: stats.gmean(x.dropna()),
                                        axis=1)

    ratio_col = "FC"
    mask = df[f"gmean_{2}"] == 0
    df[ratio_col] = df[f"gmean_{1}"] / np.where(mask, 1e-10, df[f"gmean_{2}"])

    return df


def apply_multi_group_kruskal_wallis(df: pd.DataFrame,
                                     groups: List[List[str]]) -> pd.DataFrame:
    '''
    Row-wise multi-group Kruskal-Wallis test;
    groups contains sublists of columns
    Applies scipy Kruskal-Wallis test to each row
    for groups defined by columns in 'groups',
    adds the resulting pvalue to a new column
    and returns the updated data frame
    '''
    p_values = []
    for _, row in df.iterrows():
        # Create a list of groups based on the columns in cols
        groups_values = [row[group] for group in groups]

        # Apply the Kruskal-Wallis test to the groups
        _, p_value = scipy.stats.kruskal(*groups_values, nan_policy='omit')
        p_values.append(p_value)

    df['pvalue'] = p_values
    return df


def first_column_for_column_values(df: pd.DataFrame, columns: List,
                                   values: List) -> List:
    """
    Given a dataframe df and selection columns, selects rows
    where values are equal to the those
    in the "values" List (provided in pairwise fashion).
    Returns: list of values of the first column for selected rows.
    """

    if not all(len(sublist) == len(columns) for sublist in values):
        raise ValueError(
            "Number of values in each sublist must be"
            " the same as number of columns")

    # Create a mask for each column-value pair
    first_column_values_list = []
    for vals in values:
        masks = []
        for column, value in zip(columns, vals):
            masks.append(df[column] == value)
        # Combine the masks using logical AND
        mask = reduce(lambda x, y: x & y, masks)
        first_column_values_list.append(list(df[mask].iloc[:, 0]))

    return first_column_values_list


# the name of the function is preserved but strong cleaning was performed:
#  the weird function for zero replacement no longer exists
def arg_repl_zero2value(how: str, df: pd.DataFrame) -> float:
    """
       how: string indicating how replacing zero values (e.g. "min/2")
       The result is a float
       """
    how = how.lower()
    err_msg = "replace_zero_with argument is not correctly formatted"
    if how.startswith("min"):
        if how == "min":
            denominator = int(1)  # n is the denominator, default is 1
        else:
            try:
                denominator = float(str(how.split("/")[1]))
            except Exception as e:
                logger.info(e)
                raise ValueError(err_msg)
        min_value = df[df > 0].min(skipna=True).min(skipna=True)
        output_value = min_value / denominator
    else:
        try:
            output_value = float(str(how))
        except Exception as e:
            logger.info(e)
            raise ValueError(err_msg)

    return output_value


def overlap_symmetric(x: np.array, y: np.array) -> float:
    a = [np.nanmin(x), np.nanmax(x)]
    b = [np.nanmin(y), np.nanmax(y)]

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    overlap = np.nanmax([a[0], b[0]]) - np.nanmin([a[1], b[1]])
    return overlap


def overlap_asymmetric(x: np.array, y: np.array) -> float:
    # x is the reference group
    overlap = np.nanmin(y) - np.nanmax(x)
    return overlap


def compute_distance_between_intervals(group1: np.array, group2: np.array,
                                       overlap_method: str) -> float:
    """
    computes the distance between intervals provided in group1 and group2
    """
    assert_literal(overlap_method, overlap_methods_types, "overlap method : ")

    if overlap_method == "symmetric":
        return overlap_symmetric(group1, group2)
    else:
        return overlap_asymmetric(group1, group2)


def df_to_dict_by_compartment(df: pd.DataFrame,
                              metadata: pd.DataFrame) -> dict:
    """
    splits df into a dictionary of dataframes, each for one compartment
    """
    output_dict = dict()
    for compartment in metadata["short_comp"].unique():
        sample_names = metadata[metadata["short_comp"] == compartment][
            "original_name"]
        compartment_df = df[list(sample_names)]
        output_dict[compartment] = compartment_df
    return output_dict


def check_dict_has_keys(d: dict, expected_keys: list) -> np.array:
    has_key = []
    for k in expected_keys:
        has_key.append(k in d.keys())
    return np.array(has_key)


def check_dict_has_known_values(d: dict, possible_values: list) -> np.array:
    known_val = []
    for v in d.values():
        known_val.append(v in possible_values)
    return np.array(known_val)


def verify_metadata_sample_not_duplicated(metadata_df: pd.DataFrame) -> None:
    '''
    checks for duplicated sample names in a metadata DataFrame and raises an
    error if any conflicts are detected.
    '''
    sample_counts = metadata_df["name_to_plot"].value_counts()
    duplicated_samples = sample_counts[sample_counts > 1].index.tolist()

    if duplicated_samples:
        txt_errors = f"-> duplicated sample names: {duplicated_samples}\n"
        raise ValueError(f"Found conflicts in your metadata:\n{txt_errors}")


def a12(lst1, lst2, rev=True):
    """
    Non-parametric hypothesis testing using Vargha and Delaney's
    A12 statistic:
    how often is x in lst1 greater than y in lst2?
    == > it gives a size effect, good to highlight
    potentially real effects <==
    """
    # credits : Macha Nikolski
    more = same = 0.0
    for x in lst1:
        for y in lst2:
            if x == y:
                same += 1
            elif rev and x > y:
                more += 1
            elif not rev and x < y:
                more += 1
    return (more + 0.5 * same) / (len(lst1) * len(lst2))


def compute_gmean_nonan(anarray: np.array) -> float:
    arr_nonzero = np.where(anarray == 0, np.finfo(float).eps, anarray)
    return stats.gmean(arr_nonzero[~np.isnan(arr_nonzero)])


def countnan_samples(df: pd.DataFrame, groups: List) -> pd.DataFrame:
    """
    Calculates the count of NaN values in each row of the DataFrame
    and for each group within the specified columns,
    and adds two new columns to the DataFrame with the counts.

    Only works if groups contains two sublists of column names
    """
    assert len(groups) == 2
    df["count_nan_samples_group1"] = df[groups[0]].isnull().sum(axis=1)
    df["count_nan_samples_group2"] = df[groups[1]].isnull().sum(axis=1)

    return df


def dynamic_xposition_ylabeltext(plotwidth) -> float:
    position_float = plotwidth * 0.00145
    if position_float < 0.01:
        position_float = 0.01
    return position_float


def flatten(xs):
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x


def compute_ranksums_allH0(vInterest: np.array, vBaseline: np.array):
    """
    The Wilcoxon rank-sum test tests the null hypothesis that two sets of
     measurements are drawn from the same distribution.
    ‘two-sided’: one of the distributions (underlying x or y) is
        stochastically greater than the other.
    ‘less’: the distribution underlying x is stochastically less
        than the distribution underlying y.
    ‘greater’: the distribution underlying x is stochastically
        greater than the distribution underlying y.
    """
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.ranksums(vInterest, vBaseline, alternative="less")
    sta2, p2 = scipy.stats.ranksums(vInterest, vBaseline,
                                    alternative="greater")
    sta3, p3 = scipy.stats.ranksums(vInterest, vBaseline,
                                    alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1])  # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def compute_wilcoxon_allH0(vInterest: np.array, vBaseline: np.array):
    #  Wilcoxon signed-rank test
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.wilcoxon(vInterest, vBaseline, alternative="less")
    sta2, p2 = scipy.stats.wilcoxon(vInterest, vBaseline,
                                    alternative="greater")
    sta3, p3 = scipy.stats.wilcoxon(vInterest, vBaseline,
                                    alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1])  # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def compute_brunnermunzel_allH0(vInterest: np.array, vBaseline: np.array):
    vInterest = vInterest[~np.isnan(vInterest)]
    vBaseline = vBaseline[~np.isnan(vBaseline)]
    sta, p = scipy.stats.brunnermunzel(vInterest, vBaseline,
                                       alternative="less")
    sta2, p2 = scipy.stats.brunnermunzel(vInterest, vBaseline,
                                         alternative="greater")
    sta3, p3 = scipy.stats.brunnermunzel(vInterest, vBaseline,
                                         alternative="two-sided")

    # best (smaller pvalue) among all tailed tests
    pretups = [(sta, p), (sta2, p2), (sta3, p3)]
    tups = []
    for t in pretups:  # make list of tuples with no-nan pvalues
        if not np.isnan(t[1]):
            tups.append(t)

    if len(tups) == 0:  # if all pvalues are nan assign two sided result
        tups = [(sta3, p3)]

    stap_tup = min(tups, key=lambda x: x[1])  # nan already excluded
    stat_result = stap_tup[0]
    pval_result = stap_tup[1]

    return stat_result, pval_result


def absolute_geommean_diff(b_values: np.array, a_values: np.array):
    m_b = compute_gmean_nonan(b_values)
    m_a = compute_gmean_nonan(a_values)
    diff_absolute = abs(m_b - m_a)
    return diff_absolute


def drop_all_nan_metabolites_on_comp_frames(frames_dict: Dict,
                                            metadata: pd.DataFrame) -> Dict:
    """ metabolites must be in rows """
    compartments = metadata["short_comp"].unique().tolist()
    for dataset in frames_dict.keys():
        for compartment in compartments:
            tmp = frames_dict[dataset][compartment]
            tmp = tmp.dropna(how="all", subset=tmp.columns.difference(["ID"]),
                             axis=0)
            frames_dict[dataset][compartment] = tmp
    return frames_dict


def set_samples_names(frames_dict: Dict, metadata: pd.DataFrame) -> Dict:
    """
    Given a dictionary where each dataset has been split in compartment
    dataframes,
    goes through them and renames all the columns (sample names) to those that
    we want to see on the plot; excludes the ID column from this
    """
    for dataset, compartments_dict in frames_dict.items():
        for compartment, df in compartments_dict.items():
            original_names = metadata[metadata["short_comp"] == compartment][
                "original_name"]
            new_names = metadata[metadata["short_comp"] == compartment][
                "name_to_plot"]
            renamed_columns = {old: new for old, new in
                               zip(original_names, new_names)
                               if old != "ID"}
            renamed_df = df.rename(columns=renamed_columns)
            frames_dict[dataset][compartment] = renamed_df

    return frames_dict
