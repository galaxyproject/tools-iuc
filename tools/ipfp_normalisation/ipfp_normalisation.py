#!/usr/bin/env python
"""
IPFP Normalisation
"""
import argparse
import sys

import numpy as np


def throw_error(msg, exit_code=1):
    sys.stderr.write(msg)
    sys.exit(exit_code)


def ipfp(data, precision=1e-5, maxIterations=50):
    """
    Return the normalized version of the input data (matrix) as an ndarray
    :param data:				np.ndArray
    :param precision:			float		combined allowed deviation (residual error) of col and row means from TARGET (=1)
    :param maxIterations:		int			maximum amount of iterations (1x row and 1x col per iteration)
    :return normalizedData:		np.ndArray	normalized data
    """
    try:
        assert isinstance(data, np.ndarray) and data.dtype in ['float64', 'int64']
        assert precision > 0
        assert isinstance(maxIterations, int) and maxIterations > 0
    except AssertionError:
        throw_error("Invalid input parameters. Please check that the input data consists of floats or integers, precision > 0 and maxIterations is a positive integer.")
    # replace zeros with nan
    if (data < 0).any():
        throw_error("Negative values detected, only use positive values.")

    zeros = (data == 0)
    if zeros.any():
        print("Zero values detected; replacing with NA.")
        data = data.astype(float)
        data[zeros] = np.nan

    # initialize variables
    Nrows, Ncols = data.shape
    convergenceTrail = np.asarray([np.nan] * (2 * maxIterations))
    convergence = np.inf
    normalized_data = data
    TARGET = 1

    i = 0  # number of current iteration
    # without reshaping the ndarrays, they have shape (x,) (no second value) and the procedure fails.
    # main loop; iterates until convergence is reached (i.e., L1-norm below variable <h>) or the maximum number of
    # iteration cycles is surpassed.
    while convergence > precision and i < maxIterations:
        # fit the rows
        Ri = TARGET * np.asarray(1 / np.nanmean(normalized_data, 1)).reshape(Nrows,)
        normalized_data = (normalized_data.T * Ri).T

        # calculate deviation from column marginals; row deviation is zero at even indices. (index start = 0)
        convergenceTrail[2 * i] = Nrows * 0.5 * np.nansum(np.abs(np.nanmean(normalized_data, 0) - TARGET))

        # fit the columns
        Si = TARGET * np.asarray(1 / np.nanmean(normalized_data, 0)).reshape(Ncols,)
        normalized_data *= Si
        # calculate deviation from row marginals; column deviation is zero at odd indices. (index start = 0)
        convergenceTrail[2 * i + 1] = Ncols * 0.5 * np.nansum(np.abs(np.nanmean(normalized_data, 1) - TARGET))

        convergence = convergenceTrail[2 * i + 1]
        i += 1

    if i == maxIterations:
        throw_error(f"Max number of IPFP iterations ({maxIterations}) reached. Attained precision: {convergence}.")

    return normalized_data


def main():
    parser = argparse.ArgumentParser(description="IPFP Normalisation")
    parser.add_argument('-i', '--input', help="Input file", required=True, metavar="FILE")
    parser.add_argument('-p', '--precision', help="Precision", default=1e-5, type=float)
    parser.add_argument('-m', '--maxIterations', help="Max iterations", default=50, type=int)
    parser.add_argument('-s', '--skipHeaders', help="Skip headers, skips the first n lines", default=0, type=int)

    args = parser.parse_args()

    try:
        data = np.genfromtxt(args.input, skip_header=args.skipHeaders, filling_values=np.nan, delimiter='\t')
        normalized_data = ipfp(data, args.precision, args.maxIterations)
        np.savetxt("output.tsv", normalized_data, delimiter='\t')

    except Exception as e:
        throw_error(str(e))


if __name__ == "__main__":
    main()
