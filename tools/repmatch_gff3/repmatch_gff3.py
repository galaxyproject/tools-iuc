# repmatch.py
#
# Replicate matching - matches paired peaks from two or more replicates
#
# Input: one or more gff files (matched_peak output from cwpair2, each a list of paired peaks from a replicate
#
# Output: list of matched groups and list of unmatched peaks
# Files: statistics_table.tabular (file to replicate ID), matched_paired_peaks.tabular, detail.tabular, unmatched_peaks.tabular

import argparse

import repmatch_gff3_util

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='inputs', action='append', nargs=2, help="Input datasets")
    parser.add_argument('--method', dest='method', default='closest', help='Method of finding match')
    parser.add_argument('--distance', dest='distance', type=int, default=50, help='Maximum distance between peaks in different replicates to allow merging')
    parser.add_argument('--step', dest='step', type=int, default=0, help='Step size of distance for each iteration')
    parser.add_argument('--replicates', dest='replicates', type=int, default=2, help='Minimum number of replicates that must be matched for merging to occur')
    parser.add_argument('--low_limit', dest='low_limit', type=int, default=-1000, help='Lower limit for c-w distance filter')
    parser.add_argument('--up_limit', dest='up_limit', type=int, default=1000, help='Upper limit for c-w distance filter')
    parser.add_argument('--output_files', dest='output_files', default='all', help='Restrict output dataset collections.')
    parser.add_argument('--output_matched_peaks', dest='output_matched_peaks', help='Matched groups in gff format')
    parser.add_argument('--output_unmatched_peaks', dest='output_unmatched_peaks', default=None, help='Unmatched paired peaks in tabular format')
    parser.add_argument('--output_detail', dest='output_detail', default=None, help='Details in tabular format')
    parser.add_argument('--output_statistics_table', dest='output_statistics_table', default=None, help='Keys in tabular format')
    parser.add_argument('--output_statistics_histogram', dest='output_statistics_histogram', default=None, help='Histogram')

    args = parser.parse_args()

    dataset_paths = []
    hids = []
    for (dataset_path, hid) in args.inputs:
        dataset_paths.append(dataset_path)
        hids.append(hid)
    repmatch_gff3_util.process_files(dataset_paths,
                                     hids,
                                     args.method,
                                     args.distance,
                                     args.step,
                                     args.replicates,
                                     args.up_limit,
                                     args.low_limit,
                                     args.output_files,
                                     args.output_matched_peaks,
                                     args.output_unmatched_peaks,
                                     args.output_detail,
                                     args.output_statistics_table,
                                     args.output_statistics_histogram)
