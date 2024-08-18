import argparse
import sys
from collections import Counter

import scipy


def linkage_as_newick(linkage, tip_names):
    newick_parts = tip_names[::]
    within_cluster_distances = [0] * len(tip_names)
    for step in linkage:
        n1 = int(step[0])
        n2 = int(step[1])
        d = float(step[2])
        d1 = d - within_cluster_distances[n1]
        d2 = d - within_cluster_distances[n2]
        id1 = newick_parts[n1]
        id2 = newick_parts[n2]
        part = f'({id1}:{d1 / 2},{id2}:{d2 / 2})'
        within_cluster_distances.append(d)
        newick_parts.append(part)
    return newick_parts[-1].format(*newick_parts) + ';'


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'infile',
        help='Distance matrix input file'
    )
    parser.add_argument(
        'out_prefix',
        help="Output prefix"
    )
    parser.add_argument
    parser.add_argument(
        '-m', '--method', default="average",
        choices=[
            "single",
            "complete",
            "average",
            "weighted",
            "centroid",
            "median",
            "ward"
        ],
        help="Clustering method to use"
    )
    missing_names = parser.add_mutually_exclusive_group()
    missing_names.add_argument(
        "--nc", "--no-colnames", action="store_true",
        help="Indicate that the distance matrix input does not feature column names"
    )
    missing_names.add_argument(
        "--nr", "--no-rownames", action="store_true",
        help="Indicate that the distance matrix input does not feature row names"
    )
    cut_mode = parser.add_mutually_exclusive_group()
    cut_mode.add_argument(
        "-n", "--n-clusters", nargs="*", type=int
    )
    cut_mode.add_argument(
        "--height", nargs="*", type=float
    )
    parser.add_argument("-s", "--min-cluster-size", type=int, default=2)
    args = parser.parse_args()

    # read from input and check that
    # we have been passed a symmetric distance matrix
    with open(args.infile) as i:
        col_count = None
        row_count = 0
        matrix = []
        if args.nc:
            col_names = col_count = None
        else:
            while True:
                # skip leading empty lines
                line = next(i).strip("\n\r")
                if line:
                    break
            if args.nr:
                col_names = line.split("\t")
            else:
                # first column is for row names, rest are column names
                col_names = line.split("\t")[1:]
            col_count = len(col_names)
            if not col_count:
                sys.exit(
                    'No data columns found. '
                    'By default, this tool expects tabular input with column names on the first line '
                    'and a row name in the first column of each row followed by data columns. '
                    'Use --no-colnames or --no-rownames to modify the expected format.'
                )
        for line in i:
            if not line.strip():
                # skip empty lines
                continue
            row_count += 1
            if col_count is not None and row_count > col_count:
                sys.exit(
                    'This tool expects a symmetric distance matrix with an equal number of rows and columns, '
                    'but got more rows than columns.'
                )
            if args.nr:
                row_name = None
                row_data = line.strip("\n\r").split("\t")
            else:
                row_name, *row_data = line.strip("\n\r").split("\t")
            if col_count is None:
                col_count = len(row_data)
                col_names = [None] * col_count
            col_name = col_names[row_count - 1]
            if not row_name and col_name:
                # tolerate omitted row names, use col name instead
                row_name = col_name
            elif row_name and not col_name:
                # likewise for column names
                # plus update list of col names with row name
                col_name = col_names[row_count - 1] = row_name
            elif not row_name and not col_name:
                sys.exit(
                    'Each sample in the distance matrix must have its name specified via a row name, a column name, or both, '
                    f'but found no name for sample number {row_count}'
                )
            if row_name != col_name:
                sys.exit(
                    'This tool expects a symmetric distance matrix with identical names for rows and columns, '
                    f'but got "{col_name}" in column {row_count} and "{row_name}" on row {row_count}.'
                )
            if len(row_data) != col_count:
                sys.exit(
                    'This tool expects a symmetric distance matrix with the same number of columns on each row, '
                    f'but row {row_count} ("{row_name}") has {len(row_data)} columns instead of {col_count}.'
                )
            try:
                matrix.append([float(x) for x in row_data])
            except ValueError as e:
                if args.nr:
                    sys.exit(str(e) + f' on row {row_count}')
                else:
                    sys.exit(str(e) + f' on row {row_count} ("{row_name}")')
    if row_count < col_count:
        sys.exit(
            'This tool expects a symmetric distance matrix with an equal number of rows and columns, '
            'but got more columns than rows.'
        )

    # turn the distance matrix into "condensed" vector form
    # this gives us further checks and raises ValueErrors if:
    # - the values on the diagonal aren't zero
    # - the upper and lower triangle of the matrix aren't identical
    D = scipy.spatial.distance.squareform(matrix)

    # perform the requested clustering and retrieve the result as a linkage object
    linkage = scipy.cluster.hierarchy.linkage(D, args.method)

    with open(args.out_prefix + '.tree.newick', 'w') as o:
        o.write(linkage_as_newick(linkage, col_names))

    # cut the tree as specified and report sample to cluster assignments
    if args.n_clusters or args.height:
        if args.n_clusters:
            cut_values = args.n_clusters
            colname_template = "cluster_id_n{}"
        else:
            cut_values = args.height
            colname_template = "cluster_id_h{}"
        header_cols = ["sample"] + [
            colname_template.format(x) for x in cut_values
        ]
        cut_result = scipy.cluster.hierarchy.cut_tree(
            linkage,
            args.n_clusters,
            args.height
        )

        # Go through the cut results once to determine cluster sizes

        # In the final report, the ids of clusters with fewer members than
        # args.min_cluster_size will be masked with "-".
        # The remaining cluster ids will be renumbered to start fom 1.
        # This has to be done for each clustering resulting from the
        # user-specified cut_values.
        cluster_member_counts = [Counter() for _ in cut_values]
        effective_cluster_ids = [{} for _ in cut_values]
        for cluster_ids in cut_result:
            for cl_count, cl_id, eff_id in zip(cluster_member_counts, cluster_ids, effective_cluster_ids):
                cl_count[cl_id] += 1
        for counter, eff_ids in zip(cluster_member_counts, effective_cluster_ids):
            eff_id = 1
            for item, count in counter.items():
                # Since Python 3.7, Counter objects (like dicts) preserve
                # insertion order so we can be sure that in the mapping
                # constructed below, clusters will get renumbered in
                # the order they will be reported later.
                if count >= args.min_cluster_size:
                    eff_ids[item] = str(eff_id)
                    eff_id += 1
                else:
                    eff_ids[item] = "-"

        # build and write the cluster assignment report
        # with remapped cluster ids
        cluster_assignments = []
        for name, cluster_ids in zip(col_names, cut_result):
            cluster_assignments.append(
                [name]
                + [
                    eff_ids[c]
                    for c, eff_ids in zip(cluster_ids, effective_cluster_ids)
                ]
            )
        with open(args.out_prefix + '.cluster_assignments.tsv', 'w') as o:
            print("\t".join(header_cols), file=o)
            for ass in cluster_assignments:
                print("\t".join(ass), file=o)
