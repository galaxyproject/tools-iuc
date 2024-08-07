import argparse
import sys

import scipy


def linkage_as_newick(linkage, tip_names):
    newick_parts = tip_names[::]
    within_cluster_distances = [0]*len(tip_names)
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
        choices = [
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
    cut_mode = parser.add_mutually_exclusive_group()
    cut_mode.add_argument(
        "-n", "--n-clusters", nargs="*", type=int
    )
    cut_mode.add_argument(
        "--height", nargs="*", type=float
    )
    args = parser.parse_args()

    # TO DO:
    # - parse outputs to generate

    # read from input and check that
    # we have been passed a symmetric distance matrix
    with open(args.infile) as i:
        col_names = next(i).rstrip("\n\r").split("\t")[1:]
        col_count = len(col_names)
        if not col_count:
            sys.exit(
                'No data columns found. '
                'This tool expects tabular input with column names on the first line '
                'and a row name in the first column of each row followed by data columns.'
            )
        row_count = 0
        matrix = []
        for line in i:
            if not line.strip():
                # skip empty lines
                continue
            row_count += 1
            if row_count > col_count:
                sys.exit(
                    'This tool expects a symmetric distance matrix with an equal number of rows and columns, '
                    'but got more rows than columns.'
                )
            row_name, *row_data = line.strip(" \n\r").split("\t")
            col_name = col_names[row_count - 1]
            if not row_name:
                # tolerate omitted row names, use col name instead
                row_name = col_name
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

    if args.n_clusters or args.height:
        cluster_assignments = []
        for name, cluster_ids in zip(
            col_names,
            scipy.cluster.hierarchy.cut_tree(
                linkage,
                args.n_clusters,
                args.height
            )
        ):
            cluster_assignments.append(
                [name]
                + [str(c + 1) for c in cluster_ids]
            )
        with open(args.out_prefix + '.cluster_assignments.tsv', 'w') as o:
            for ass in cluster_assignments:
                print("\t".join(ass), file=o)
