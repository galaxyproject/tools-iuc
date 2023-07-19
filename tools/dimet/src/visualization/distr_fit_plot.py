import logging
import operator
import os
from functools import reduce
from typing import List

from constants import (
    assert_literal,
    availtest_methods_type,
    data_files_keys_type,
)

from data import Dataset

import helpers

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

from omegaconf import DictConfig

import pandas as pd

from processing import fit_statistical_distribution
from processing.differential_analysis import \
    select_rows_with_sufficient_non_nan_values

import scipy.stats as stats

logger = logging.getLogger(__name__)


def make_pdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function"""

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if \
        arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if \
        arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf


def plot_best_fit(data, dist_str, pdf, out_file) -> None:
    plt.figure(figsize=(12, 8))
    plt.hist(data, bins="auto", density=True, alpha=0.5, label="Data")
    plt.plot(pdf, lw=2, label="PDF")
    plt.legend(loc="upper right", shadow=True, fontsize="x-large")
    plt.title("Best fit distribution \n" + dist_str)
    plt.xlabel("z-score")
    plt.ylabel("frequency")
    # plt.set_title(u'Best fit distribution \n' + dist_str)
    plt.savefig(out_file)
    logger.info(f"saved plot to {out_file}")


def find_best_distribution_to_plot(df: pd.DataFrame, out_file):
    """
    Find the best distribution among all the scipy.stats distributions
    and return it together with its parameters
    """
    logger.info("Fitting a distribution")
    dist = np.around(np.array((df["zscore"]).astype(float)), 5)

    best_dist, best_dist_name, best_fit_params = get_best_fit_to_plot(
        dist, out_file)

    logger.info(f"Best fit is {best_dist_name} with {best_fit_params}")
    args_param = dict(e.split("=") for e in best_fit_params.split(", "))
    for k, v in args_param.items():
        args_param[k] = float(v)

    best_distribution = getattr(stats, best_dist_name)
    q_val = best_dist.ppf(0.95, **args_param)
    logger.info(f"And the q value is {q_val}")
    return best_distribution, args_param


def get_best_fit_to_plot(input_array, out_file):
    matplotlib.rcParams["figure.figsize"] = (16.0, 12.0)
    matplotlib.style.use("ggplot")
    """Return the best fit distribution to data and its parameters"""

    # Load data
    data = pd.Series(input_array)

    # Find best fit distribution
    best_fit_name, best_fit_params = \
        fit_statistical_distribution.best_fit_distribution(data, 200)

    best_dist = getattr(stats, best_fit_name)

    # Make probability density function (PDF) with best params
    pdf = make_pdf(best_dist, best_fit_params)

    # parameters
    param_names = (best_dist.shapes + ", loc, scale").split(", ") if \
        best_dist.shapes else ["loc", "scale"]
    param_str = ", ".join(["{}={:0.2f}".format(k, v) for k, v
                           in zip(param_names, best_fit_params)])

    # Display
    dist_str = '{} ({})'.format(best_fit_name, param_str)
    plot_best_fit(data, dist_str, pdf, out_file)

    return best_dist, best_fit_name, param_str


def run_dist_fit_plot_pairwise(
    df: pd.DataFrame, dataset: Dataset, cfg: DictConfig,
    comparison: List[str], test: availtest_methods_type, out_file_path: str
) -> None:
    """
    Runs a pairwise comparison when distribution fitting test
    for plotting only
    """
    assert test == "disfit"
    conditions_list = helpers.first_column_for_column_values(
        df=dataset.metadata_df, columns=cfg.analysis.method.grouping,
        values=comparison
    )
    # flatten the list of lists and select the subset of column names
    # present in the sub dataframe
    columns = [i for i in reduce(operator.concat, conditions_list)
               if i in df.columns]
    this_comparison = [list(filter(lambda x: x in columns, sublist))
                       for sublist in conditions_list]
    df4c = df[columns].copy()
    df4c = df4c[(df4c.T != 0).any()]  # delete rows being zero everywhere
    df4c = df4c.dropna(axis=0, how="all")
    df4c = helpers.row_wise_nanstd_reduction(df4c)
    df4c = helpers.countnan_samples(df4c, this_comparison)

    df4c = helpers.calculate_gmean(df4c, this_comparison)
    df_good, df_bad = select_rows_with_sufficient_non_nan_values(
        df4c, groups=this_comparison)

    df_good = fit_statistical_distribution.compute_z_score(df_good, "FC")

    find_best_distribution_to_plot(df_good, out_file_path)


def run_distr_fit_plot(
        file_name: data_files_keys_type, dataset: Dataset,
        cfg: DictConfig, test: availtest_methods_type,
        out_plot_dir: str, mode: str) -> None:
    """
    Differential comparison is performed
    Attention: we replace zero values using the provided method
    distribution fitting plots are saved
    """
    assert_literal(test, availtest_methods_type, "Available test")
    assert_literal(file_name, data_files_keys_type, "file name")

    impute_value = cfg.analysis.method.impute_values[file_name]
    for compartment, compartmentalized_df in \
            dataset.compartmentalized_dfs[file_name].items():
        df = compartmentalized_df
        df = df[(df.T != 0).any()]
        val_instead_zero = helpers.arg_repl_zero2value(impute_value,
                                                       df)
        df = df.replace(to_replace=0, value=val_instead_zero)
        if mode == "pairwise":
            for comparison in cfg.analysis.comparisons:
                comp = "-".join(map(lambda x: "-".join(x), comparison))
                file_basename = dataset.get_file_for_label(file_name)
                file_basename += f"--{compartment}-{comp}-{test}.pdf"
                out_file_path = os.path.join(out_plot_dir, file_basename)
                run_dist_fit_plot_pairwise(df, dataset, cfg, comparison,
                                           test, out_file_path)
        # elif mode == "time_course":
        #     pass  # not implemented to date, evaluate if worthed
        logger.info(f"saved plots to {out_plot_dir}")
