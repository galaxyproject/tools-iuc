from scipy.special import gammaln as log_gamma, logsumexp as log_sum_exp, psi

import numba
import numpy as np


def fit_annealed(
    log_p_data,
    priors,
    var_params,
    annealing_power=1.0,
    convergence_threshold=1e-6,
    max_iters=int(1e4),
    num_annealing_steps=10,
    print_freq=100,
):
    if num_annealing_steps == 1:
        annealing_ladder = [1.0]

    else:
        annealing_ladder = np.linspace(0, 1.0, num_annealing_steps) ** annealing_power

    for t in annealing_ladder:
        print("Setting annealing factor to : {}".format(t))
        print()

        log_p_data_annealed = t * log_p_data

        if t == 1.0:
            convergence_threshold_t = convergence_threshold

        else:
            convergence_threshold_t = convergence_threshold * 1e-2

        elbo_trace = fit(
            log_p_data_annealed,
            priors,
            var_params,
            convergence_threshold=convergence_threshold_t,
            max_iters=max_iters,
            print_freq=print_freq,
        )

    return elbo_trace


def fit(
    log_p_data,
    priors,
    var_params,
    convergence_threshold=1e-6,
    max_iters=int(1e4),
    print_freq=100,
):
    elbo_trace = [compute_elbo(log_p_data, priors, var_params)]

    for i in range(max_iters):
        if i % print_freq == 0:
            num_clusters = len(set(var_params.z.argmax(axis=1)))
            print("Iteration: {}".format(i))
            print("ELBO: {}".format(elbo_trace[-1]))
            print("Number of clusters used: {}".format(num_clusters))
            print()

        update_z(log_p_data, var_params)

        update_pi(priors, var_params)

        update_theta(log_p_data, priors, var_params)

        elbo_trace.append(compute_elbo(log_p_data, priors, var_params))

        diff = (elbo_trace[-1] - elbo_trace[-2]) / np.abs(elbo_trace[-1])

        if diff < convergence_threshold:
            break

    return elbo_trace


def get_priors(num_clusters, num_grid_points):
    return Priors(
        np.ones(num_clusters), (1 / num_grid_points) * np.ones(num_grid_points)
    )


def get_variational_params(num_clusters, num_data_points, num_dims, num_grid_points):
    var_params = VariationalParameters(
        np.random.dirichlet(np.ones(num_clusters)),
        np.random.gamma(1, 1, size=(num_clusters, num_dims, num_grid_points)),
        np.random.dirichlet(np.ones(num_clusters), size=num_data_points),
    )

    var_params.theta = (
        var_params.theta / np.sum(var_params.theta, axis=2)[:, :, np.newaxis]
    )

    return var_params


class Priors(object):
    def __init__(self, pi, theta):
        self.pi = pi

        self.theta = theta


class VariationalParameters(object):
    def __init__(self, pi, theta, z):
        self.pi = pi

        self.theta = theta

        self.z = z


def compute_elbo(log_p_data, priors, var_params):
    return compute_e_log_p(log_p_data, priors, var_params) - compute_e_log_q(var_params)


def compute_e_log_p(log_p_data, priors, var_params):
    log_p = 0

    log_p += log_gamma(np.sum(priors.pi)) - np.sum(log_gamma(priors.pi))

    log_p += np.sum(
        (priors.pi + np.sum(var_params.z, axis=0) - 1)
        * (psi(var_params.pi) - psi(np.sum(var_params.pi)))
    )

    log_p += np.sum(var_params.theta * np.log(priors.theta)[np.newaxis, np.newaxis, :])

    log_p += np.sum(
        var_params.z * compute_log_p_data_theta(log_p_data, var_params.theta)
    )

    return log_p


def compute_e_log_q(var_params):
    log_p = 0

    log_p += log_gamma(np.sum(var_params.pi)) - np.sum(log_gamma(var_params.pi))

    log_p += np.sum(
        (var_params.pi - 1) * (psi(var_params.pi) - psi(np.sum(var_params.pi)))
    )

    log_p += np.sum(var_params.theta * np.log(var_params.theta + 1e-6))

    log_p += np.sum(var_params.z * np.log(var_params.z + 1e-6))

    return log_p


def update_pi(priors, var_params):
    var_params.pi = priors.pi + np.sum(var_params.z, axis=0)


def update_z(log_p_data, var_params):
    var_params.z = compute_log_p_data_theta(log_p_data, var_params.theta)

    var_params.z += (psi(var_params.pi) - psi(np.sum(var_params.pi)))[np.newaxis, :]

    var_params.z = var_params.z - log_sum_exp(var_params.z, axis=1)[:, np.newaxis]

    var_params.z = np.exp(var_params.z)


def update_theta(log_p_data, priors, var_params):
    var_params.theta = np.log(
        priors.theta[np.newaxis, np.newaxis, :]
    ) + compute_log_p_data_z(log_p_data, var_params.z)

    var_params.theta = (
        var_params.theta - log_sum_exp(var_params.theta, axis=2)[:, :, np.newaxis]
    )

    var_params.theta = np.exp(var_params.theta)


@numba.njit(parallel=True)
def compute_log_p_data_z(log_p_data, z):
    """Equivalent to np.sum(var_params.z[:, :, np.newaxis, np.newaxis] * log_p_data[:, np.newaxis, :, :], axis=0)"""
    N, D, G = log_p_data.shape

    K = z.shape[1]

    result = np.zeros((K, D, G))

    for k in numba.prange(K):
        for d in range(D):
            for g in range(G):
                for n in range(N):
                    result[k, d, g] += log_p_data[n, d, g] * z[n, k]

    return result


@numba.njit(parallel=True)
def compute_log_p_data_theta(log_p_data, theta):
    """Equivalent to np.sum(var_params.theta[np.newaxis, :, :, :] * log_p_data[:, np.newaxis, :, :], axis=(2, 3))"""
    N, D, G = log_p_data.shape

    K = theta.shape[0]

    result = np.zeros((N, K))

    for n in numba.prange(N):
        for k in range(K):
            for d in range(D):
                for g in range(G):
                    result[n, k] += log_p_data[n, d, g] * theta[k, d, g]

    return result
