import sys
import pandas as pd


def main():
    inp = sys.argv[1]
    df = pd.read_csv(inp, sep='\t', dtype=str)

    lab_idx = int(sys.argv[2])
    count_idx = int(sys.argv[3])

    ncols = len(df.columns)
    if not (0 <= lab_idx < ncols) or not (0 <= count_idx < ncols):
        print(f"ERROR: column index out of range. File has {ncols} columns.",
              file=sys.stderr)
        sys.exit(2)

    xy_column = df.columns[lab_idx]
    count_column = df.columns[count_idx]

    # Always normalize to sets (handles NaN/empty safely)
    df[count_column] = (
        df[count_column]
        .fillna("")
        .astype(str)
        .str.split(",")
        .map(lambda xs: {x.strip() for x in xs if x and x.strip()})
    )

    x = count_column + "_x"
    y = count_column + "_y"
    cross = df[[xy_column, count_column]].merge(df[[xy_column, count_column]],
                                                how="cross")

    cross["intersection"] = (
        cross.apply(lambda row: row[x].intersection(row[y]), axis=1)
        .map(lambda s: ",".join(sorted(s)) if s else None)
    )

    c_mat = pd.DataFrame(index=df[xy_column], columns=df[xy_column])

    x = xy_column + "_x"
    y = xy_column + "_y"
    mat = (cross.pivot(index=x, columns=y, values="intersection")
           .rename_axis(None, axis=1).rename_axis(None))
    mat = mat.reindex(df[xy_column]).reindex(df[xy_column], axis=1)

    for i in mat:
        c_mat[i] = (mat[i]
                    .fillna("")
                    .astype(str)
                    .str.split(",")
                    .apply(
            lambda xs: len([x.strip() for x in xs if x and x.strip()]))
                    )

    mat.to_csv(sys.argv[4], sep="\t")
    c_mat.to_csv(sys.argv[5], sep="\t")


if __name__ == "__main__":
    main()
