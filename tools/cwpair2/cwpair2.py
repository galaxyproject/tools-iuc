"""
cwpair2.py

Takes a list of called peaks on both strands and produces a list of matched pairs and a list
of unmatched orphans using a specified method for finding matched pairs.  Methods for finding
matched pairs are mode, closest, largest or all, where the analysis is run for each method

Input: list of one or more gff format files

Output: files produced for each input/mode combination:
MP (matched_pair), D (details), O (orphans), P (frequency preview plot), F (frequency final plot),
C (statistics graph), statistics.tabular
"""

import argparse
import csv

import cwpair2_util

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='inputs', action='append', nargs=2, help="Input datasets")
    parser.add_argument('--method', dest='method', default='mode', help='Method of finding match.')
    parser.add_argument('--up_distance', dest='up_distance', type=int, default=50, help='Distance upstream from a pair.')
    parser.add_argument('--down_distance', dest='down_distance', type=int, default=100, help='Distance downstream of a pair.')
    parser.add_argument('--binsize', dest='binsize', type=int, default=1, help='Width of bins for plots and mode.')
    parser.add_argument('--threshold_format', dest='threshold_format', help='Percentage to filter the 95th percentile.')
    parser.add_argument('--relative_threshold', dest='relative_threshold', type=float, default=0.0, help='Percentage to filter the 95th percentile.')
    parser.add_argument('--absolute_threshold', dest='absolute_threshold', type=float, default=0.0, help='Absolute value to filter.')
    parser.add_argument('--output_files', dest='output_files', default='matched_pair', help='Restrict output dataset collections.')
    parser.add_argument('--statistics_output', dest='statistics_output', help='Statistics output file.')
    args = parser.parse_args()

    cwpair2_util.create_directories()

    statistics = []
    if args.absolute_threshold > 0:
        threshold = args.absolute_threshold
    elif args.relative_threshold > 0:
        threshold = args.relative_threshold / 100.0
    else:
        threshold = 0
    for (dataset_path, hid) in args.inputs:
        stats = cwpair2_util.process_file(dataset_path,
                                          hid,
                                          args.method,
                                          threshold,
                                          args.up_distance,
                                          args.down_distance,
                                          args.binsize,
                                          args.output_files)
        statistics.extend(stats)
    # Accumulate statistics.
    by_file = {}
    for stats in statistics:
        # Skip "None" statistics from failed files
        if not stats:
            continue
        path = stats['stats_path']
        if path not in by_file:
            by_file[path] = []
        by_file[path].append(stats)
    # Write tabular statistics file.
    keys = ['fname', 'final_mode', 'preview_mode', 'perc95', 'paired', 'orphans']
    statistics_out = csv.writer(open(args.statistics_output, 'wt'), delimiter='\t')
    statistics_out.writerow(keys)
    for file_path, statistics in by_file.items():
        for stats in statistics:
            statistics_out.writerow([stats[key] for key in keys])
