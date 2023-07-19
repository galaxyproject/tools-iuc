#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
import logging
import warnings

import numpy as np

import pandas as pd

import scipy.stats as stats


logger = logging.getLogger(__name__)


def compute_z_score(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """
    Add one column with z-score of ratio already computed
    """
    df = df.assign(zscore=stats.zscore(df[column_name]))  # was "ratio"
    return df


def find_best_distribution(df: pd.DataFrame):
    """
    Find the best distribution among all the scipy.stats distributions
    and return it together with its parameters
    """
    logger.info("Fitting a distribution")
    dist = np.around(np.array((df["zscore"]).astype(float)), 5)

    best_dist, best_dist_name, best_fit_params = get_best_fit(dist)

    logger.info(f"Best fit is {best_dist_name} with {best_fit_params}")
    args_param = dict(e.split("=") for e in best_fit_params.split(", "))
    for k, v in args_param.items():
        args_param[k] = float(v)

    best_distribution = getattr(stats, best_dist_name)
    q_val = best_dist.ppf(0.95, **args_param)
    logger.info(f"And the q value is {q_val}")
    return best_distribution, args_param


def get_best_fit(input_array):
    """Return the best fit distribution to data and its parameters"""

    # Load data
    data = pd.Series(input_array)

    # Find best fit distribution
    best_fit_name, best_fit_params = best_fit_distribution(data, 200)

    best_dist = getattr(stats, best_fit_name)

    # parameters
    param_names = (best_dist.shapes + ", loc, scale").split(", ") if \
        best_dist.shapes else ["loc", "scale"]
    param_str = ", ".join(["{}={:0.2f}".format(k, v)
                           for k, v in zip(param_names, best_fit_params)])

    return best_dist, best_fit_name, param_str


def best_fit_distribution(data, bins=200):
    # TODO: try and catch exception and warnings
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)

    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Distributions to check
    # TODO: Get distribution list not hardcoded ?
    DISTRIBUTIONS = [
        stats.alpha,
        stats.anglit,
        stats.arcsine,
        stats.beta,
        stats.betaprime,
        stats.bradford,
        stats.burr,
        stats.cauchy,
        stats.chi,
        stats.chi2,
        stats.cosine,
        stats.dgamma,
        stats.dweibull,
        stats.erlang,
        stats.expon,
        stats.exponnorm,
        stats.exponweib,
        stats.exponpow,
        stats.f,
        stats.fatiguelife,
        stats.foldcauchy,
        stats.foldnorm,
        stats.genlogistic,
        stats.genpareto,
        stats.gennorm,
        stats.genexpon,
        stats.genextreme,
        stats.gausshyper,
        stats.gamma,
        stats.gengamma,
        stats.genhalflogistic,
        stats.gilbrat,
        stats.gompertz,
        stats.gumbel_r,
        stats.gumbel_l,
        stats.halfcauchy,
        stats.halflogistic,
        stats.halfnorm,
        stats.halfgennorm,
        stats.hypsecant,
        stats.invgamma,
        stats.invgauss,
        stats.invweibull,
        stats.johnsonsb,
        stats.johnsonsu,
        stats.ksone,
        stats.kstwobign,
        stats.laplace,
        stats.levy,
        stats.levy_l,
        stats.fisk,
        stats.logistic,
        stats.loggamma,
        stats.loglaplace,
        stats.lognorm,
        stats.lomax,
        stats.maxwell,
        stats.mielke,
        stats.nakagami,
        stats.ncx2,
        stats.ncf,
        stats.nct,
        stats.norm,
        stats.pareto,
        stats.pearson3,
        stats.powerlaw,
        stats.powerlognorm,
        stats.powernorm,
        stats.rdist,
        stats.reciprocal,
        stats.rayleigh,
        stats.rice,
        stats.recipinvgauss,
        stats.semicircular,
        stats.t,
        stats.triang,
        stats.truncexpon,
        stats.truncnorm,
        stats.tukeylambda,
        stats.uniform,
        stats.vonmises,
        stats.vonmises_line,
        stats.wald,
        stats.weibull_min,
        stats.weibull_max,
        stats.wrapcauchy,
    ]

    # Best holders
    best_distribution = stats.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in DISTRIBUTIONS:
        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))
                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse

        except Exception:
            pass

    return best_distribution.name, best_params
