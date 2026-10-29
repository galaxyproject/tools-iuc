import json
import optparse
import os
import re
import tempfile
import unicodedata

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from b2bTools import MultipleSeq

# Columns contributed by each predictor. The tabular export of b2bTools uses a
# fixed schema that also holds the PSPer columns, which stay empty unless that
# predictor runs, so only the columns of the selected predictors are kept.
PREDICTOR_COLUMNS = {
    "dynamine": ["backbone", "sidechain", "helix", "sheet", "coil", "ppII"],
    "disomine": ["disoMine"],
    "efoldmine": ["earlyFolding"],
    "agmata": ["agmata"],
}
# DynaMine is always executed, and DisoMine and AgMata are built on the
# output of EFoldMine, so b2bTools runs these whether or not they were asked
# for. Their columns are reported rather than dropped.
ALWAYS_EXECUTED = ["dynamine"]
REQUIRED_PREDICTORS = {
    "disomine": ["efoldmine"],
    "agmata": ["efoldmine"],
}
PREDICTOR_ORDER = ["dynamine", "disomine", "efoldmine", "agmata"]
INDEX_COLUMNS = ["sequence_id", "residue_index", "residue"]
# Biophysical regions and value ranges, as used by the single sequence tool.
THRESHOLDS = {
    "backbone": {
        "membrane spanning": [1.0, 1.5],
        "rigid": [0.8, 1.0],
        "context-dependent": [0.69, 0.8],
        "flexible": [-1.0, 0.69],
    },
    "earlyFolding": {
        "early folds": [0.169, 2.0],
        "late folds": [-1.0, 0.169],
    },
    "disoMine": {"ordered": [-1.0, 0.5], "disordered": [0.5, 2.0]},
}
ORDERED_REGIONS = {
    "backbone": [
        "flexible",
        "context-dependent",
        "rigid",
        "membrane spanning",
    ],
    "earlyFolding": ["late folds", "early folds"],
    "disoMine": ["ordered", "disordered"],
}
REGION_COLORS = ["yellow", "orange", "pink", "red"]
VALUE_RANGES = {
    "backbone": [-0.2, 1.2],
    "sidechain": [-0.2, 1.2],
    "ppII": [-0.2, 1.2],
    "earlyFolding": [-0.2, 1.2],
    "disoMine": [-0.2, 1.2],
    "agmata": [-0.2, 1.2],
    "helix": [-1.0, 1.0],
    "sheet": [-1.0, 1.0],
    "coil": [-1.0, 1.0],
}
DISTRIBUTION_COLUMNS = {
    "median": "median",
    "firstQuartile": "first_quartile",
    "thirdQuartile": "third_quartile",
    "bottomOutlier": "bottom_outlier",
    "topOutlier": "top_outlier",
}


def slugify(value):
    """
    From https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode'. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    value = (
        unicodedata.normalize("NFKD", value)
        .encode("ascii", "ignore")
        .decode("ascii")
    )
    value = re.sub(r"[^\w\s-]", "", value.lower())
    return re.sub(r"[-\s]+", "-", value).strip("-_")


def rounder_function(value):
    return round(float(value), 3)


def selected_tools(options):
    """List the predictors b2bTools will actually execute."""
    requested = set(ALWAYS_EXECUTED)
    for tool in ("disomine", "efoldmine", "agmata"):
        if getattr(options, tool):
            requested.add(tool)
            requested.update(REQUIRED_PREDICTORS.get(tool, []))
    return [tool for tool in PREDICTOR_ORDER if tool in requested]


def build_predictions_dataframe(msa, tools):
    """Read the aligned per-residue predictions as a data frame.

    The alignment mapping (including the gap positions, which carry no
    predicted value) is delegated to b2bTools itself instead of being derived
    again here.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        tabular_path = os.path.join(tmp_dir, "predictions.tsv")
        msa.get_all_predictions_tabular(tabular_path, sep="\t")
        dataframe = pd.read_csv(tabular_path, sep="\t")
    prediction_columns = sorted(
        column for tool in tools for column in PREDICTOR_COLUMNS[tool]
    )
    dataframe = dataframe[INDEX_COLUMNS + prediction_columns]
    return dataframe.round(decimals=3)


def build_distribution_dataframe(distributions):
    """Flatten the per-alignment-position statistics into a tidy table.

    One row per alignment position and predicted feature keeps the set of
    columns independent of the selected predictors, which makes the result
    straightforward to filter and group with the other Galaxy tools.
    """
    rows = []
    for prediction in sorted(distributions.keys()):
        statistics = distributions[prediction]
        positions = len(statistics["median"])
        for position in range(positions):
            row = {
                "alignment_position": position,
                "prediction": prediction,
            }
            for source_name, column_name in DISTRIBUTION_COLUMNS.items():
                row[column_name] = statistics[source_name][position]
            rows.append(row)
    dataframe = pd.DataFrame(
        rows,
        columns=[
            "alignment_position",
            "prediction",
            *DISTRIBUTION_COLUMNS.values(),
        ],
    )
    return dataframe.round(decimals=3)


def write_json(predictions, output_filepath):
    payload = {
        "sequences": predictions["sequences"],
        "proteins": {
            sequence_key: {
                key: value
                for key, value in sequence_predictions.items()
                if isinstance(value, list)
            }
            for sequence_key, sequence_predictions in predictions[
                "proteins"
            ].items()
        },
    }
    rounded_payload = json.loads(
        json.dumps(payload), parse_float=rounder_function
    )
    with open(output_filepath, "w") as json_file:
        json_file.write(json.dumps(rounded_payload, indent=2, sort_keys=True))


def split_by_sequence(dataframe, output_dir):
    for sequence_key, sequence_df in dataframe.groupby(
        "sequence_id", sort=False
    ):
        filename = os.path.join(output_dir, f"{slugify(sequence_key)}.tsv")
        sequence_df.drop(columns=["sequence_id"]).to_csv(
            filename, sep="\t", index=False
        )


def as_float_array(values):
    """Turn a list of predicted values into a plottable array.

    Alignment positions where no sequence has a value are ``None``, which
    matplotlib cannot handle, so they become ``NaN`` and leave a gap in the
    chart instead.
    """
    return np.array(
        [np.nan if value is None else float(value) for value in values],
        dtype=float,
    )


def check_min_max(values, former_min, former_max):
    """Widen a value range so that the plotted values always fit in it."""
    known_values = values[np.isfinite(values)]
    if known_values.size == 0:
        return former_min, former_max
    if known_values.max() + 0.1 > former_max:
        former_max = float(known_values.max()) + 0.1
    if known_values.min() - 0.1 < former_min:
        former_min = float(known_values.min()) - 0.1
    return former_min, former_max


def place_legend(ax, legend):
    """Keep a legend that sits outside the axes visible.

    add_artist() clips its artist to the axes, which would hide a legend
    anchored outside of them, so the clipping is switched off.
    """
    ax.add_artist(legend)
    legend.set_clip_on(False)


def highlight_regions(ax, prediction):
    """Shade the biophysical regions defined for a predicted feature.

    Only a few features have meaningful thresholds; the others are left
    without any shading.
    """
    if prediction not in ORDERED_REGIONS:
        return
    for index, region in enumerate(ORDERED_REGIONS[prediction]):
        lower, upper = THRESHOLDS[prediction][region]
        ax.axhspan(
            lower, upper, alpha=0.3, color=REGION_COLORS[index], label=region
        )
    handles, labels = ax.get_legend_handles_labels()
    handles_by_label = dict(zip(labels, handles))
    # Sorted "from up to low" so the legend matches the shading.
    regions = list(reversed(ORDERED_REGIONS[prediction]))
    place_legend(
        ax,
        ax.legend(
            [handles_by_label[region] for region in regions],
            regions,
            fancybox=True,
            shadow=True,
            loc="lower left",
            bbox_to_anchor=(1.04, 0),
        ),
    )


def plot_distribution(prediction, statistics, plot_output, highlight):
    series = {
        name: as_float_array(values) for name, values in statistics.items()
    }
    positions = range(len(series["median"]))
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    ax.set_title(f"{prediction} distribution along the alignment")
    min_value, max_value = VALUE_RANGES[prediction]
    for name in series:
        min_value, max_value = check_min_max(series[name], min_value, max_value)
    ax.fill_between(
        positions,
        series["bottomOutlier"],
        series["topOutlier"],
        alpha=0.15,
        color="tab:blue",
        label="outlier boundaries",
    )
    ax.fill_between(
        positions,
        series["firstQuartile"],
        series["thirdQuartile"],
        alpha=0.35,
        color="tab:blue",
        label="interquartile range",
    )
    ax.plot(
        positions, series["median"], color="tab:blue", label="median"
    )
    place_legend(
        ax,
        ax.legend(
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
            fancybox=True,
            shadow=True,
        ),
    )
    if highlight:
        highlight_regions(ax, prediction)
    ax.set_xlim([0, len(series["median"]) - 1])
    ax.set_ylim([min_value, max_value])
    ax.set_xlabel("alignment position")
    ax.set_ylabel("prediction values")
    ax.grid(axis="y")
    plt.savefig(
        os.path.join(plot_output, f"{slugify(prediction)}.png"),
        bbox_inches="tight",
    )
    plt.close()


def plot_aligned_values(prediction, sequences_values, plot_output, highlight):
    """Draw the predicted values of every sequence over the alignment.

    Each sequence is one line, so the same alignment position can be
    compared across sequences. Gap positions have no value, so the line of
    a sequence is interrupted wherever it has a gap.
    """
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    ax.set_title(f"{prediction} prediction by aligned position")
    min_value, max_value = VALUE_RANGES[prediction]
    alignment_length = 0
    for sequence_key in sorted(sequences_values.keys()):
        values = as_float_array(sequences_values[sequence_key])
        min_value, max_value = check_min_max(values, min_value, max_value)
        ax.plot(range(len(values)), values, label=sequence_key)
        alignment_length = max(alignment_length, len(values))
    place_legend(
        ax,
        ax.legend(
            bbox_to_anchor=(1.04, 1),
            loc="upper left",
            fancybox=True,
            shadow=True,
        ),
    )
    if highlight:
        highlight_regions(ax, prediction)
    ax.set_xlim([0, alignment_length - 1])
    ax.set_ylim([min_value, max_value])
    ax.set_xlabel("alignment position")
    ax.set_ylabel("prediction values")
    ax.grid(axis="y")
    plt.savefig(
        os.path.join(plot_output, f"aligned_{slugify(prediction)}.png"),
        bbox_inches="tight",
    )
    plt.close()


def main(options):
    tools = selected_tools(options)
    msa = MultipleSeq()
    msa.from_aligned_file(options.input_msa, tools=tools)

    predictions_df = build_predictions_dataframe(msa, tools)
    predictions_df.to_csv(options.tabular_output, sep="\t", index=False)
    if options.split_output:
        split_by_sequence(predictions_df, options.split_output)

    distributions = msa.get_all_predictions_msa_distrib()["results"]
    distribution_df = build_distribution_dataframe(distributions)
    distribution_df.to_csv(
        options.distribution_output, sep="\t", index=False
    )

    predictions = msa.get_all_predictions_msa()
    write_json(predictions, options.json_output)

    prediction_keys = sorted(
        column for tool in tools for column in PREDICTOR_COLUMNS[tool]
    )
    if options.plot_distribution:
        for prediction in prediction_keys:
            plot_distribution(
                prediction,
                distributions[prediction],
                options.plot_output,
                options.highlight,
            )
    if options.plot_all:
        for prediction in prediction_keys:
            sequences_values = {
                sequence_key: sequence_predictions[prediction]
                for sequence_key, sequence_predictions in predictions[
                    "proteins"
                ].items()
            }
            plot_aligned_values(
                prediction,
                sequences_values,
                options.plot_output,
                options.highlight,
            )


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option(
        "--disomine",
        action="store_true"
    )
    parser.add_option(
        "--efoldmine",
        action="store_true"
    )
    parser.add_option(
        "--agmata",
        action="store_true"
    )
    parser.add_option(
        "--file",
        dest="input_msa",
        type="string"
    )
    parser.add_option(
        "--tabular",
        dest="tabular_output",
        type="string"
    )
    parser.add_option(
        "--distribution",
        dest="distribution_output",
        type="string"
    )
    parser.add_option(
        "--split-output",
        dest="split_output",
        type="string"
    )
    parser.add_option(
        "--json",
        dest="json_output",
        type="string"
    )
    parser.add_option(
        "--plot-output",
        dest="plot_output",
        type="string"
    )
    parser.add_option(
        "--plot_distribution",
        action="store_true"
    )
    parser.add_option(
        "--plot_all",
        action="store_true"
    )
    parser.add_option(
        "--highlight",
        action="store_true"
    )
    try:
        options, args = parser.parse_args()
        if not options.input_msa:
            parser.error('Input file not given (--file)')
        if not options.tabular_output:
            parser.error('Tabular output file not given (--tabular)')
        if not options.distribution_output:
            parser.error('Distribution output file not given (--distribution)')
        if not options.json_output:
            parser.error('Json output file not given (--json)')
        if (options.plot_distribution or options.plot_all) and not options.plot_output:
            parser.error('Plot output directory not given (--plot-output)')
        main(options)
    except optparse.OptionError as exc:
        raise RuntimeError(f"Invalid arguments: {args}") from exc
