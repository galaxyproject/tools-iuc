import math
import re
import sys

import bokeh.palettes
import numpy as np
import pandas as pd
from bokeh.io import save
from bokeh.models import ColorBar, ColumnDataSource, LinearColorMapper
from bokeh.plotting import figure
from bokeh.transform import transform


def main():
    # Expected args:
    # 1) counts matrix (TSV)
    # 2) names matrix (TSV)
    # 3) threshold (int)
    # 4) output.tsv HTML
    # 5) x-axis label (optional string)
    # 6) y-axis label (optional string)
    if len(sys.argv) < 5:
        print(
            "ERROR: missing arguments. Expected: <count_mat> <name_mat> "
            "<threshold> <out.html> [x_label] [y_label]",
            file=sys.stderr,
        )
        sys.exit(2)

    # Read matrices. Force dtype=str to preserve labels exactly as strings.
    c_raw = pd.read_csv(sys.argv[1], sep="\t", index_col=0, dtype=str)
    n_mat = pd.read_csv(sys.argv[2], sep="\t", index_col=0, dtype=str)

    # Convert count matrix to numeric (float), keep NaN if conversion fails.
    c_mat = c_raw.apply(pd.to_numeric, errors="coerce").astype(float)

    # Basic sanity check: same row/column labels in both matrices.
    if not c_mat.index.equals(n_mat.index) or not c_mat.columns.equals(n_mat.columns):
        print(
            "ERROR: input matrices must have the same row/column labels.",
            file=sys.stderr,
        )
        sys.exit(2)

    # Threshold
    threshold = int(sys.argv[3])

    # Set diagonal to NaN to avoid bias in display counting
    np.fill_diagonal(c_mat.values, float("nan"))

    # Filtering phase (use numeric count matrix, not tuples/strings)
    if threshold != 0:
        rem = []
        for col in c_mat.columns:
            col_max = np.nanmax(c_mat[col].values)
            if np.isnan(col_max) or col_max < threshold:
                rem.append(col)

        if rem:
            c_mat.drop(rem, inplace=True, axis=1)
            c_mat.drop(rem, inplace=True, axis=0)
            n_mat.drop(rem, inplace=True, axis=1)
            n_mat.drop(rem, inplace=True, axis=0)

    # Build long-form dataframe for Bokeh:
    # - idx: row label
    # - col: column label
    # - value: numeric count
    # - names: cleaned string list
    df_counts = c_mat.stack(dropna=False).rename("value").reset_index()
    df_counts.columns = ["idx", "col", "value"]

    df_names = n_mat.stack(dropna=False).rename("names").reset_index()
    df_names.columns = ["idx", "col", "names"]

    df = pd.merge(df_counts, df_names, on=["idx", "col"], how="left")

    # Clean names (remove brackets/quotes)
    df["names"] = (
        df["names"].fillna("").astype(str).map(lambda s: re.sub(r"[\[\]']", "", s))
    )

    # Ensure numeric values
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    # Use a continuous mapper to avoid categorical palette length issues
    palette = bokeh.palettes.inferno(256)
    if df["value"].notna().any():
        vmin = float(np.nanmin(df["value"].values))
        vmax = float(np.nanmax(df["value"].values))
        if vmin == vmax:
            # Avoid a degenerate mapper if all values are identical
            vmax = vmin + 1.0
    else:
        vmin, vmax = 0.0, 1.0

    mapper = LinearColorMapper(palette=palette, low=vmin, high=vmax, nan_color="gray")

    # Axis labels (optional)
    x_label = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != "" else ""
    y_label = sys.argv[6] if len(sys.argv) > 6 and sys.argv[6] != "" else ""

    # Define a figure
    p = figure(
        width=800,
        height=800,
        title="Heatmap",
        x_axis_label=x_label,
        y_axis_label=y_label,
        x_range=list(df["idx"].drop_duplicates()),
        y_range=list(df["col"].drop_duplicates()[::-1]),
        tooltips=[
            ("Column", "@col"),
            ("Index", "@idx"),
            ("# of common elements", "@value"),
            ("Common elements", "@names"),
        ],
        x_axis_location="above",
        output_backend="webgl",
    )

    # IMPORTANT: Bokeh expects radians, not degrees
    p.xaxis.major_label_orientation = math.pi / 4

    # Create rectangles for heatmap
    source = ColumnDataSource(df)
    p.rect(
        x="idx",
        y="col",
        width=1,
        height=1,
        source=source,
        line_color=None,
        fill_color=transform("value", mapper),
    )

    # Add color bar
    color_bar = ColorBar(
        color_mapper=mapper,
        location=(0, 0),
        label_standoff=6,
        border_line_color=None,
    )
    p.add_layout(color_bar, "right")

    # Save HTML file
    save(p, filename=sys.argv[4])


if __name__ == "__main__":
    main()
