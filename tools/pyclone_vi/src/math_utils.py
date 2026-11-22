import numba
import numpy as np


# Special functions


@numba.njit
def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)


@numba.njit
def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(n - x) - log_factorial(x)


@numba.njit
def log_factorial(x):
    return log_gamma(x + 1)


@numba.njit
def log_gamma(x):
    return np.math.lgamma(x)


@numba.njit
def log_normalize(x):
    return x - log_sum_exp(x)


@numba.njit
def log_sum_exp(log_X):
    max_exp = np.max(log_X)

    if np.isinf(max_exp):
        return max_exp

    total = 0

    for x in log_X:
        total += np.exp(x - max_exp)

    return np.log(total) + max_exp


# Probability densities


@numba.njit
def log_beta_binomial_pdf(n, x, a, b):
    return log_binomial_coefficient(n, x) + log_beta(a + x, b + n - x) - log_beta(a, b)


@numba.njit
def log_binomial_pdf(n, x, p):
    return log_binomial_coefficient(n, x) + x * np.log(p) + (n - x) * np.log1p(-p)
