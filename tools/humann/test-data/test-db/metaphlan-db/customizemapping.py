#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script to generate a extract a custom mapping file from input mapping file.
# Mostly used for a reduced-size demo data generation.


import argparse
from pathlib import Path


if __name__ == '__main__':
    # Read command line
    parser = argparse.ArgumentParser(description='Customize HUMAnN utility mapping')
    parser.add_argument('--in_mapping', help="Path to mapping file to reduce")
    parser.add_argument('--features', help="Path to tabular file with features to keep in first column")
    parser.add_argument('--elements', help="Path to tabular file with elements to keep in other columns")
    parser.add_argument('--out_mapping', help="Path to reduced mapping file")
    args = parser.parse_args()

    in_mapping_fp = Path(args.in_mapping)
    feature_fp = Path(args.features)
    element_fp = Path(args.elements)
    out_mapping_fp = Path(args.out_mapping)

    # extract features to keep
    features = set()
    with open(feature_fp, 'r') as feature_f:
        for line in feature_f.readlines():
            features.add(line.split("\t")[0])
    print(features)

    # extract elements to keep
    elements = set()
    with open(element_fp, 'r') as element_f:
        for line in element_f.readlines():
            elements.add(line.split("\t")[0])
    print(elements)

    # write mapping for features to keep while keeping only elements
    with open(in_mapping_fp, 'r') as in_mapping_f:
        with open(out_mapping_fp, 'w') as out_mapping_f:
            for line in in_mapping_f.readlines():
                l_split = line.split("\t")
                feat = l_split[0]
                if feat in features:
                    to_write = [feat]
                    for e in l_split[1:]:
                        if e in elements:
                            to_write.append(e)
                    out_mapping_f.write("%s\n" % '\t'.join(to_write))
