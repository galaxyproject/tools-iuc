import argparse
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import plotly.subplots as sp
import statsmodels.api as sm  # to build a LOWESS model
from scipy.stats import gaussian_kde


# subplot titles
def make_subplot_titles(sample_names: List[str]) -> List[str]:
    """Generates subplot titles for the MA plot.

    Args:
        sample_names (list): List of sample names.

    Returns:
        list: List of subplot titles.
    """
    subplot_titles = []
    num_samples = len(sample_names)
    for i in range(num_samples):
        for j in range(num_samples):
            if i == j:
                subplot_titles.append(f"{sample_names[i]}")
            else:
                subplot_titles.append(f"{sample_names[i]} vs. {sample_names[j]}")
    return subplot_titles


def densities(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Calculates the density of points for a scatter plot.

    Args:
        x (array-like): X-axis values.
        y (array-like): Y-axis values.

    Returns:
        array: Density values for the points.
    """
    values = np.vstack([x, y])
    return gaussian_kde(values)(values)


def movingaverage(data: np.ndarray, window_width: int) -> np.ndarray:
    """Calculates the moving average of the data.

    Args:
        data (array-like): Input data.
        window_width (int): Width of the moving window.

    Returns:
        array: Moving average values.
    """
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    return ma_vec


def update_max(current: float, values: np.ndarray) -> float:
    """Updates the maximum value.

    Args:
        current (float): Current maximum value.
        values (array-like): Array of values to compare.

    Returns:
        float: Updated maximum value.
    """
    return max(current, np.max(values))


def get_indices(
    num_samples: int, num_cols: int, plot_num: int
) -> Tuple[int, int, int, int]:
    """Calculates the indices for subplot placement.

    Args:
        num_samples (int): Number of samples.
        num_cols (int): Number of columns in the subplot grid.
        plot_num (int): Plot number.

    Returns:
        tuple: Indices for subplot placement (i, j, col, row).
    """
    i = plot_num // num_samples
    j = plot_num % num_samples
    col = plot_num % num_cols + 1
    row = plot_num // num_cols + 1
    return i, j, col, row


def create_subplot_data(
    frac: float,
    it: int,
    num_bins: int,
    window_width: int,
    samples: pd.DataFrame,
    i: int,
    j: int,
) -> Dict:
    """Creates data for a single subplot.

    Args:
        frac (float): LOESS smoothing parameter.
        it (int): Number of iterations for LOESS smoothing.
        num_bins (int): Number of bins for histogram.
        window_width (int): Window width for moving average.
        samples (DataFrame): DataFrame containing sample data.
        i (int): Index of the first sample.
        j (int): Index of the second sample.

    Returns:
        dict: Data for the subplot.
    """
    subplot_data = {}
    subplot_data["mean"] = np.log(samples.iloc[:, [i, j]].mean(axis=1))
    if i == j:
        counts, bins = np.histogram(subplot_data["mean"], bins=num_bins)
        subplot_data["bins"] = bins
        subplot_data["counts"] = counts
        subplot_data["counts_smoothed"] = movingaverage(counts, window_width)
        subplot_data["max_counts"] = np.max(counts)
    else:
        subplot_data["log_fold_change"] = np.log2(
            samples.iloc[:, i] / samples.iloc[:, j]
        )
        subplot_data["max_log_fold_change"] = np.max(subplot_data["log_fold_change"])
        subplot_data["densities"] = densities(
            subplot_data["mean"], subplot_data["log_fold_change"]
        )
        subplot_data["regression"] = sm.nonparametric.lowess(
            subplot_data["log_fold_change"], subplot_data["mean"], frac=frac, it=it
        )
    return subplot_data


def create_plot_data(
    frac: float,
    it: int,
    num_bins: int,
    window_width: int,
    samples: pd.DataFrame,
    num_samples: int,
    num_plots: int,
    num_cols: int,
) -> List[Dict]:
    """Creates data for all subplots.

    Args:
        frac (float): LOESS smoothing parameter.
        it (int): Number of iterations for LOESS smoothing.
        num_bins (int): Number of bins for histogram.
        window_width (int): Window width for moving average.
        samples (DataFrame): DataFrame containing sample data.
        num_samples (int): Number of samples.
        num_plots (int): Number of plots.
        num_cols (int): Number of columns in the subplot grid.

    Returns:
        list: List of data for each subplot.
    """
    plots_data = []
    for plot_num in range(num_plots):
        i, j, _, _ = get_indices(num_samples, num_cols, plot_num)
        subplot_data = create_subplot_data(
            frac, it, num_bins, window_width, samples, i, j
        )
        plots_data.append(subplot_data)
    return plots_data


def ma_plots_plotly(
    num_rows: int,
    num_cols: int,
    num_plots: int,
    plots_data: List[Dict],
    sample_names: List[str],
    size: int,
    ylim_hist: float,
    ylim_ma: float,
    features: np.ndarray,
) -> go.Figure:
    """Generates MA plots using Plotly.

    Args:
        num_rows (int): Number of rows in the subplot grid.
        num_cols (int): Number of columns in the subplot grid.
        num_plots (int): Number of plots.
        plots_data (list): List of data for each subplot.
        sample_names (list): List of sample names.
        size (int): Size of the plot.
        ylim_hist (float): Y-axis limit for histograms.
        ylim_ma (float): Y-axis limit for MA plots.
        features (array-like): Feature names.

    Returns:
        Figure: Plotly figure object.
    """
    fig = sp.make_subplots(
        rows=num_rows,
        cols=num_cols,
        shared_xaxes="all",
        subplot_titles=make_subplot_titles(sample_names),
    )

    for plot_num in range(num_plots):
        i, j, col, row = get_indices(len(sample_names), num_cols, plot_num)
        subplot_data = plots_data[plot_num]

        mean = subplot_data["mean"]

        if i == j:
            # Plot histogram on the diagonal
            hist_bar = go.Bar(
                x=subplot_data["bins"],
                y=subplot_data["counts"],
            )
            fig.add_trace(hist_bar, row=row, col=col)

            hist_line = go.Scatter(
                x=subplot_data["bins"],
                y=subplot_data["counts_smoothed"],
                marker=dict(
                    color="red",
                ),
            )
            fig.add_trace(hist_line, row=row, col=col)
            fig.update_yaxes(
                title_text="Counts",
                range=[0, ylim_hist],
                matches="y1",
                showticklabels=True,
                row=row,
                col=col,
            )
        else:
            log_fold_change = subplot_data["log_fold_change"]
            scatter = go.Scatter(
                x=mean,
                y=log_fold_change,
                mode="markers",
                marker=dict(
                    color=subplot_data["densities"], symbol="circle", colorscale="jet"
                ),
                name=f"{sample_names[i]} vs {sample_names[j]}",
                text=features,
                hovertemplate="<b>%{text}</b><br>Log Mean: %{x}<br>Log2 Fold Change: %{y}<extra></extra>",
            )
            fig.add_trace(scatter, row=row, col=col)

            regression = subplot_data["regression"]
            line = go.Scatter(
                x=regression[:, 0],
                y=regression[:, 1],
                mode="lines",
                line=dict(color="red"),
                name=f"LOWESS {sample_names[i]} vs. {sample_names[j]}",
            )
            fig.add_trace(line, row=row, col=col)

            fig.update_yaxes(
                title_text="Log2 Fold Change",
                range=[-ylim_ma, ylim_ma],
                matches="y2",
                showticklabels=True,
                row=row,
                col=col,
            )
        fig.update_xaxes(
            title_text="Log Mean Intensity", showticklabels=True, row=row, col=col
        )

    # Update layout for the entire figure
    fig.update_layout(
        height=size * num_rows,
        width=size * num_cols,
        showlegend=False,
        template="simple_white",  # Apply the 'plotly_white' template
    )
    return fig


def ma_plots_matplotlib(
    num_rows: int,
    num_cols: int,
    num_plots: int,
    pots_data: List[Dict],
    sample_names: List[str],
    size: int,
    ylim_hist: float,
    ylim_ma: float,
    window_width: int,
) -> plt.Figure:
    """Generates MA plots using Matplotlib.

    Args:
        num_rows (int): Number of rows in the subplot grid.
        num_cols (int): Number of columns in the subplot grid.
        num_plots (int): Number of plots.
        pots_data (list): List of data for each subplot.
        sample_names (list): List of sample names.
        size (int): Size of the plot.
        ylim_hist (float): Y-axis limit for histograms.
        ylim_ma (float): Y-axis limit for MA plots.
        window_width (int): Window width for moving average.

    Returns:
        Figure: Matplotlib figure object.
    """
    subplot_titles = make_subplot_titles(sample_names)
    fig, axes = plt.subplots(
        num_rows,
        num_cols,
        figsize=(size * num_cols / 100, size * num_rows / 100),
        dpi=300,
        sharex="all",
    )
    axes = axes.flatten()

    for plot_num in range(num_plots):
        i, j, _, _ = get_indices(len(sample_names), num_cols, plot_num)
        subplot_data = pots_data[plot_num]

        mean = subplot_data["mean"]

        ax = axes[plot_num]

        if i == j:
            # Plot histogram on the diagonal
            ax.bar(
                subplot_data["bins"][:-1],
                subplot_data["counts"],
                width=np.diff(subplot_data["bins"]),
                edgecolor="black",
                align="edge",
            )

            # Plot moving average line
            ax.plot(
                subplot_data["bins"][window_width // 2: -window_width // 2],
                subplot_data["counts_smoothed"],
                color="red",
            )

            ax.set_ylabel("Counts")
            ax.set_ylim(0, ylim_hist)
        else:
            # Scatter plot
            ax.scatter(
                mean,
                subplot_data["log_fold_change"],
                c=subplot_data["densities"],
                cmap="jet",
                edgecolor="black",
                label=f"{sample_names[i]} vs {sample_names[j]}",
            )

            # Regression line
            regression = subplot_data["regression"]
            ax.plot(
                regression[:, 0],
                regression[:, 1],
                color="red",
                label=f"LOWESS {sample_names[i]} vs. {sample_names[j]}",
            )

            ax.set_ylabel("Log2 Fold Change")
            ax.set_ylim(-ylim_ma, ylim_ma)

        ax.set_xlabel("Log Mean Intensity")
        ax.tick_params(labelbottom=True)  # Force showing x-tick labels
        ax.set_title(subplot_titles[plot_num])  # Add subplot title

    # Adjust layout
    plt.tight_layout()
    return fig


def main():
    """Main function to generate MA plots."""
    parser = argparse.ArgumentParser(description="Generate MA plots.")
    parser.add_argument("--file_path", type=str, help="Path to the input CSV file")
    parser.add_argument("--file_extension", type=str, help="File extension")
    parser.add_argument(
        "--frac", type=float, default=4 / 5, help="LOESS smoothing parameter"
    )
    parser.add_argument(
        "--it", type=int, default=5, help="Number of iterations for LOESS smoothing"
    )
    parser.add_argument(
        "--num_bins", type=int, default=100, help="Number of bins for histogram"
    )
    parser.add_argument(
        "--window_width", type=int, default=5, help="Window width for moving average"
    )
    parser.add_argument("--size", type=int, default=500, help="Size of the plot")
    parser.add_argument(
        "--scale", type=int, default=3, help="Scale factor for the plot"
    )
    parser.add_argument(
        "--y_scale_factor", type=float, default=1.1, help="Y-axis scale factor"
    )
    parser.add_argument(
        "--max_num_cols",
        type=int,
        default=100,
        help="Maximum number of columns in the plot",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Generate interactive plot using Plotly",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="pdf",
        choices=["pdf", "png", "html"],
        help="Output format for the plot",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="ma_plot",
        help="Output file name without extension",
    )

    args = parser.parse_args()

    # Load the data
    file_extension = args.file_extension.lower()
    if file_extension == "csv":
        data = pd.read_csv(args.file_path)
    elif file_extension in ["txt", "tsv", "tabular"]:
        data = pd.read_csv(args.file_path, sep="\t")
    elif file_extension == "parquet":
        data = pd.read_parquet(args.file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_extension}")

    features = data.iloc[:, 0]  # Assuming the first column is the feature names
    samples = data.iloc[:, 1:]  # and the rest are samples

    # Create a subplot figure
    num_samples = samples.shape[1]
    sample_names = samples.columns
    num_plots = num_samples**2
    num_cols = min(num_samples, args.max_num_cols)
    num_rows = int(np.ceil(num_plots / num_cols))

    plots_data = create_plot_data(
        args.frac,
        args.it,
        args.num_bins,
        args.window_width,
        samples,
        num_samples,
        num_plots,
        num_cols,
    )

    count_max = np.max([x.get("max_counts", 0) for x in plots_data])
    log_fold_change_max = np.max([x.get("max_log_fold_change", 0) for x in plots_data])

    ylim_hist = count_max * args.y_scale_factor
    ylim_ma = log_fold_change_max * args.y_scale_factor

    if args.interactive:
        fig = ma_plots_plotly(
            num_rows,
            num_cols,
            num_plots,
            plots_data,
            sample_names,
            args.size,
            ylim_hist,
            ylim_ma,
            features,
        )
        fig.show()
        if args.output_format == "html":
            fig.write_html(f"{args.output_file}")
        else:
            pio.write_image(
                fig,
                f"{args.output_file}",
                format=args.output_format,
                width=args.size * num_cols,
                height=args.size * num_rows,
                scale=args.scale,
            )
    else:
        fig = ma_plots_matplotlib(
            num_rows,
            num_cols,
            num_plots,
            plots_data,
            sample_names,
            args.size,
            ylim_hist,
            ylim_ma,
            args.window_width,
        )
        plt.show()
        fig.savefig(f"{args.output_file}", format=args.output_format, dpi=300)
    return 0


if __name__ == "__main__":
    main()
