#!/usr/bin/env python3.6
"""Calculate scores for Heinz.

This script transform a p-value into a score:
    1. Use alpha and lambda to calculate a threshold P-value.
    2. Calculate a score based on each P-value by alpha and the threshold.

For more details, please refer to the paper doi:10.1093/bioinformatics/btn161

Input:
    P-values from DESeq2 result: first column: names, second column P-values
Output:
    Scores, which will be used as the input of Heinz.
    First column: names, second column: scores.

Python 3.5+ is required.
"""
# Implemented by: Chao (Cico) Zhang 
# Homepage: https://Hi-IT.org
# Date: 14 Mar 2017
# Last modified: 23 May 2018

import argparse
import sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Transform a P-value into a '
                                 'score which can be used as the input of '
                                 'Heinz')
parser.add_argument('-n', '--node', required=True, dest='nodes',
                    metavar='nodes_pvalue.txt', type=str,
                    help='Input file of nodes with P-values')
parser.add_argument('-f', '--fdr', required=True, dest='fdr',
                    metavar='0.007', type=float, help='Choose a value of FDR')
parser.add_argument('-m', '--model', required=False, dest='param_file',
                    metavar='param.txt', type=str,
                    help='A txt file contains model params as input')
parser.add_argument('-a', '--alpha', required=False, dest='alpha',
                    metavar='0.234', type=float, default=0.5,
                    help='Single parameter alpha as input if txt input is '
                    'not provided')
parser.add_argument('-l', '--lambda', required=False, dest='lam',
                    metavar='0.345', type=float, default=0.5,
                    help='Single parameter lambda as input if txt input is '
                    'not provided')
parser.add_argument('-o', '--output', required=True, dest='output',
                    metavar='scores.txt', type=str,
                    help='The output file to store the calculated scores')
args = parser.parse_args()

# Check if the parameters are complete
if args.output is None:
    sys.exit('Output file is not designated.')

if args.nodes is None:
    sys.exit('Nodes with p-values must be provided.')

if args.fdr is None:
    sys.exit('FDR must be provided')

if args.fdr >= 1 or args.fdr <= 0:
    sys.exit('FDR must greater than 0 and smaller than 1')

# run heinz-print according to the input type
if args.param_file is not None:  # if BUM output is provided
    with open(args.param_file) as p:
        params = p.readlines()
        lam = float(params[0])  # Maybe this is a bug
        alpha = float(params[1])  # Maybe this is a bug
# if BUM output is not provided
elif args.alpha is not None and args.lam is not None:
    lam = args.lam
    alpha = args.alpha
else:  # The input is not complete
    sys.exit('The parameters of the model are incomplete.')

# Calculate the threshold P-value
pie = lam + (1-lam) * alpha
p_threshold = np.power((pie - lam * args.fdr)/(args.fdr - lam * args.fdr),
                       1/(alpha - 1))
print(p_threshold)
# Calculate the scores
input_pvalues = pd.read_csv(args.nodes, sep='\t', names=['node', 'pvalue'])
input_pvalues.loc[:, 'score'] = input_pvalues.pvalue.apply(lambda x:
                                                           (alpha - 1) *
                                                           (np.log(x) -
                                                            np.log(
                                                                p_threshold)))
# print(input_pvalues.loc[:, ['node', 'score']])
input_pvalues.loc[:, ['node', 'score']].to_csv(args.output, sep='\t',
                                               index=False, header=False)