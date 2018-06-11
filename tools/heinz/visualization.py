#!/usr/bin/env python3.6
"""Visualise the output of Heinz.

This script is used to visualize the output of Heinz, which is in the form of
DOT language:
    1. Clear the output of Heinz, extract the DOT source code.
    2. Visualize the DOT source code and save it into file.

The function of this script is rather simple, for more advanced visualization,
please adopt other solutions mentioned in the paper
doi: 10.1093/bioinformatics/btv526

This tool is only designed for visualizing the output of Heinz tool.
"""

# Author: Cico Zhang
# Date: 2 Aug 2017
# Version: 0.2

from graphviz import Source
import argparse
import sys


def get_args():
    """Collect the inputs."""
    parser = argparse.ArgumentParser(
        description='Visualise the output of Heinz')
    parser.add_argument('-i', '--input', required=True, dest='heinz',
                        metavar='Heinz_output.txt', type=str,
                        help='Output file of Heinz as the input')
    parser.add_argument('-o', '--output', required=True, dest='output',
                        metavar='graph.pdf', type=str,
                        help='The output file that saves the visualisation')
    args = parser.parse_args()

    if args.heinz is None:
        sys.exit('Input file must be designated.')

    return args


def main():
    """Main function."""
    args = get_args()
    # Read the whole output file
    with open(args.heinz) as r:
        graph_dot = r.readlines()

    # Remove the redundant lines
    while not graph_dot[0].startswith('graph G {'):
        graph_dot.pop(0)

    src = Source(''.join(graph_dot))
    data_pdf = src.pipe('pdf')
    # Redirect the output (very important)
    with open(args.output, 'wb') as w:
        w.write(data_pdf)
    print('The visualization is saved as PDF!')
    sys.exit(0)


if __name__ == "__main__":
    main()
