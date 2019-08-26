#!/usr/bin/python

import getopt
import math
import sys

import numpy as np

"""
Program to calculate the contrast thresholds for heatmap from tagPileUp CDT
"""


def rebin(a, new_shape):
    M, N = a.shape
    m, n = new_shape
    if m >= M:
        # repeat rows in data matrix
        a = np.repeat(a, math.ceil(float(m) / M), axis=0)

    M, N = a.shape
    m, n = new_shape

    row_delete_num = M % m
    col_delete_num = N % n

    np.random.seed(seed=0)

    if row_delete_num > 0:
        # select deleted rows with equal intervals
        row_delete = np.linspace(0, M - 1, num=row_delete_num, dtype=int)
        # sort the random selected deleted row ids
        row_delete = np.sort(row_delete)
        row_delete_plus1 = row_delete[1:-1] + \
            1  # get deleted rows plus position
        # get deleted rows plus position (top +1; end -1)
        row_delete_plus1 = np.append(
            np.append(row_delete[0] + 1, row_delete_plus1), row_delete[-1] - 1)
        # put the info of deleted rows into the next rows by mean
        a[row_delete_plus1, :] = (
            a[row_delete, :] + a[row_delete_plus1, :]) / 2
        a = np.delete(a, row_delete, axis=0)  # random remove rows

    if col_delete_num > 0:
        # select deleted cols with equal intervals
        col_delete = np.linspace(0, N - 1, num=col_delete_num, dtype=int)
        # sort the random selected deleted col ids
        col_delete = np.sort(col_delete)
        col_delete_plus1 = col_delete[1:-1] + \
            1  # get deleted cols plus position
        # get deleted cols plus position (top +1; end -1)
        col_delete_plus1 = np.append(
            np.append(col_delete[0] + 1, col_delete_plus1), col_delete[-1] - 1)
        # put the info of deleted cols into the next cols by mean
        a[:, col_delete_plus1] = (
            a[:, col_delete] + a[:, col_delete_plus1]) / 2
        a = np.delete(a, col_delete, axis=1)  # random remove columns

    M, N = a.shape

    # compare the heatmap matrix
    a_compress = a.reshape((m, M / m, n, N / n)).mean(3).mean(1)
    return np.array(a_compress)


def load_Data(input_file, quantile, absolute, header, start_col, row_num, col_num, min_upper_lim):
    data = open(input_file, 'r')
    if header == 'T':
        data.readline()

    data0 = []

    for rec in data:
        tmp = [(x.strip()) for x in rec.split('\t')]
        data0.append(tmp[start_col:])
    data0 = np.array(data0, dtype=float)

    if row_num == -999:
        row_num = data0.shape[0]
    if col_num == -999:
        col_num = data0.shape[1]

    # rebin data0
    if row_num < data0.shape[0] and col_num < data0.shape[1]:
        data0 = rebin(data0, (row_num, col_num))
    elif row_num < data0.shape[0]:
        data0 = rebin(data0, (row_num, data0.shape[1]))
    elif col_num < data0.shape[1]:
        data0 = rebin(data0, (data0.shape[0], col_num))

    # Calculate contrast limits here
    rows, cols = np.nonzero(data0)
    upper_lim = np.percentile(data0[rows, cols], quantile)
    lower_lim = 0
    if absolute != -999:
        upper_lim = absolute

    # Setting an absolute threshold to a minimum,
    # in cases the 95th percentile contrast is <= user defined min_upper_lim
    if quantile != 90.0:
        print "\nQUANTILE: {}".format(quantile)
        print "Quantile calculated UPPER LIM: {}".format(upper_lim)
        print "LOWER LIM: {}".format(lower_lim)
        if upper_lim <= min_upper_lim:
            print "setting heatmap upper_threshold to min_upper_lim\n"
            upper_lim = min_upper_lim

    outfile = open('calcThreshold.txt', 'w')
    outfile.write("upper_threshold:{}\nlower_threshold:{}\nrow_num:{}\ncol_num:{}\nheader:{}\nstart_col:{}".format(
        upper_lim, lower_lim, row_num, col_num, header, start_col))
    print('heatmap_upper_threshold:' + str(upper_lim))
    print('heatmap_lower_threshold:' + str(lower_lim))
    outfile.flush()
    outfile.close()


############################################################################
# python cdt_to_heatmap.py -i test.tabular.split_line -o test.tabular.split_line.png -q 0.9 -c black -d T -s 2 -r 500 -l 300 -b test.colorsplit
############################################################################
usage = """
Usage:
This script calculates the contrast thresholds from Tag pile up heatmap data. Outputs a text file that contains the parameters for the heatmap script.

python calculateThreshold.py -i <input file> -q <quantile> -m <min upper thresold after quantile calculation> -t <absolute tag threshold> -d <header T/F> -s <start column> -r <row num after compress> -l <col num after compress>'

Example:
python calculateThreshold.py -i test.tabular.split_line -q 90 -m 5 -d T -s 2 -r 600 -l 300
"""

if __name__ == '__main__':

    # check for command line arguments
    if len(sys.argv) < 2 or not sys.argv[1].startswith("-"):
        sys.exit(usage)
    # get arguments
    try:
        optlist, alist = getopt.getopt(sys.argv[1:], 'hi:o:q:t:c:d:s:r:l:m:')
    except getopt.GetoptError:
        sys.exit(usage)

    # default quantile contrast saturation = 0.9
    quantile = 90.0
    min_upper_lim = 5.0
    # absolute contrast saturation overrides quantile
    absolute = -999

    # default figure width/height is defined by matrix size
    # if user-defined size is smaller than matrix, activate rebin function
    row_num = -999
    col_num = -999

    for opt in optlist:
        if opt[0] == "-h":
            sys.exit(usage)
        elif opt[0] == "-i":
            input_file = opt[1]
        elif opt[0] == "-q":
            quantile = float(opt[1])
        elif opt[0] == '-t':
            absolute = float(opt[1])
        elif opt[0] == "-d":
            header = opt[1]
        elif opt[0] == "-s":
            start_col = int(opt[1])
        elif opt[0] == "-r":
            row_num = int(opt[1])
        elif opt[0] == "-l":
            col_num = int(opt[1])
        elif opt[0] == "-m":
            min_upper_lim = float(opt[1])

    print "Header present:", header
    print "Start column:", start_col
    print "Row number (pixels):", row_num
    print "Col number (pixels):", col_num
    print "Min Upper Limit while using Quantile:", min_upper_lim
    if absolute != -999:
        print "Absolute tag contrast threshold:", absolute
    else:
        print "Percentile tag contrast threshold:", quantile

    if absolute == -999 and quantile <= 0:
        print "\nInvalid threshold!!!"
        sys.exit(usage)

    load_Data(input_file, quantile, absolute,
              header, start_col, row_num, col_num, min_upper_lim)
