import json
import optparse
import os.path
import re
import unicodedata

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from b2bTools import SingleSeq


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


def check_min_max(predicted_values, former_min, former_max):
    seq_max = max(predicted_values)
    seq_min = min(predicted_values)
    if (
        seq_max + 0.1 > former_max
        and not np.isnan(seq_max)
        and not np.isinf(seq_max)
    ):
        former_max = seq_max + 0.1
    if (
        seq_min - 0.1 < former_min
        and not np.isnan(seq_min)
        and not np.isinf(seq_min)
    ):
        former_min = seq_min - 0.1
    return former_min, former_max


def plot_prediction(pred_name, hlighting_regions, predicted_values, seq_name):
    thresholds_dict = {
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
    ordered_regions_dict = {
        "backbone": [
            "flexible",
            "context-dependent",
            "rigid",
            "membrane spanning",
        ],
        "earlyFolding": ["late folds", "early folds"],
        "disoMine": ["ordered", "disordered"],
    }
    colors = ["yellow", "orange", "pink", "red"]
    ranges_dict = {
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
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    ax.set_title(pred_name + " " + "prediction")
    min_value, max_value = ranges_dict[pred_name]
    if seq_name == "all":
        max_len = 0
        for seq in predicted_values.keys():
            predictions = predicted_values[seq]
            min_value, max_value = check_min_max(
                predictions, min_value, max_value
            )
            ax.plot(range(len(predictions)), predictions, label=seq)
            if len(predictions) > max_len:
                max_len = len(predictions)
            ax.set_xlim([0, max_len - 1])
    else:
        predictions = predicted_values
        min_value, max_value = check_min_max(predictions, min_value, max_value)
        ax.plot(range(len(predictions)), predictions, label=seq_name)
        ax.set_xlim([0, len(predictions) - 1])
    legend_lines = plt.legend(
        bbox_to_anchor=(1.04, 1), loc="upper left", fancybox=True, shadow=True
    )
    ax.add_artist(legend_lines)
    # Define regions
    if hlighting_regions:
        if pred_name in ordered_regions_dict.keys():
            for i, prediction in enumerate(ordered_regions_dict[pred_name]):
                lower = thresholds_dict[pred_name][prediction][0]
                upper = thresholds_dict[pred_name][prediction][1]
                color = colors[i]
                ax.axhspan(
                    lower, upper, alpha=0.3, color=color, label=prediction
                )
            included_in_regions_legend = list(
                reversed(
                    [
                        prediction
                        for prediction in ordered_regions_dict[pred_name]
                    ]
                )
            )  # to sort it "from up to low"
            # Get handles and labels
            handles, labels = plt.gca().get_legend_handles_labels()
            handles_dict = {
                label: handles[idx] for idx, label in enumerate(labels)
            }
            # Add legend for regions, if available
            region_legend = ax.legend(
                [
                    handles_dict[region]
                    for region in included_in_regions_legend
                ],
                [region for region in included_in_regions_legend],
                fancybox=True,
                shadow=True,
                loc="lower left",
                bbox_to_anchor=(1.04, 0),
            )
            ax.add_artist(region_legend)
    ax.set_ylim([min_value, max_value])
    ax.set_xlabel("residue index")
    ax.set_ylabel("prediction values")
    ax.grid(axis="y")
    plt.savefig(
        os.path.join(
            options.plot_output,
            "{0}_{1}.png".format(slugify(seq_name), pred_name),
        ),
        bbox_inches="tight",
    )
    plt.close()


def df_dict_to_dict_of_values(df_dict, predictor):
    results_dict = {}
    for seq in df_dict.keys():
        df = pd.read_csv(df_dict[seq], sep="\t")
        results_dict[seq] = df[predictor]
    return results_dict


def main(options):
    single_seq = SingleSeq(options.input_fasta)
    b2b_tools = []
    if options.dynamine:
        b2b_tools.append("dynamine")
    if options.disomine:
        b2b_tools.append("disomine")
    if options.efoldmine:
        b2b_tools.append("efoldmine")
    if options.agmata:
        b2b_tools.append("agmata")
    single_seq.predict(b2b_tools)
    predictions = single_seq.get_all_predictions()

    def rounder_function(value):
        return round(float(value), 3)

    rounded_predictions = json.loads(
        json.dumps(predictions), parse_float=rounder_function
    )
    results_json = json.dumps(rounded_predictions, indent=2, sort_keys=True)
    with open(options.json_output, "w") as f:
        f.write(results_json)
    first_sequence_key = next(iter(predictions))
    prediction_keys = predictions[first_sequence_key].keys()
    # Sort column names
    tsv_column_names = list(prediction_keys)
    tsv_column_names.remove("seq")
    tsv_column_names = ['residue', *sorted(tsv_column_names)]

    df_dictionary = {}
    for sequence_key, seq_preds in predictions.items():
        residues = seq_preds["seq"]
        residues_count = len(residues)
        sequence_df = pd.DataFrame(
            columns=prediction_keys, index=range(residues_count)
        )
        sequence_df.index.name = "residue_index"
        for predictor in prediction_keys:
            sequence_df[predictor] = seq_preds[predictor]
        sequence_df = sequence_df.rename(columns={"seq": "residue"})
        sequence_df = sequence_df.round(decimals=3)
        filename = f"{options.output}/{slugify(sequence_key)}.tsv"
        df_dictionary[sequence_key] = filename
        sequence_df.to_csv(
            filename,
            header=True,
            columns=tsv_column_names,
            sep="\t"
        )
        # Plot each individual plot (compatible with plot all)
        if options.plot:
            for predictor in prediction_keys:
                if predictor != "seq":
                    plot_prediction(
                        pred_name=predictor,
                        hlighting_regions=True,
                        predicted_values=seq_preds[predictor],
                        seq_name=sequence_key,
                    )
    # Plot all together (compatible with plot individual)
    if options.plot_all:
        for predictor in prediction_keys:
            if predictor != "seq":
                results_dictionary = df_dict_to_dict_of_values(
                    df_dictionary, predictor
                )
                plot_prediction(
                    pred_name=predictor,
                    hlighting_regions=True,
                    predicted_values=results_dictionary,
                    seq_name="all",
                )


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option(
        "--dynamine",
        action="store_true"
    )
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
        dest="input_fasta",
        type="string"
    )
    parser.add_option(
        "--output",
        dest="output",
        type="string"
    )
    parser.add_option(
        "--plot-output",
        type="string",
        dest="plot_output"
    )
    parser.add_option(
        "--json",
        dest="json_output",
        type="string"
    )
    parser.add_option(
        "--plot",
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
        if not (options.dynamine or options.disomine or options.efoldmine or options.agmata):
            parser.error('At least one predictor is required')
        if not options.input_fasta:
            parser.error('Input file not given (--file)')
        if not options.output:
            parser.error('Output directory not given (--output)')
        if (options.plot or options.plot_all) and not options.plot_output:
            parser.error('Plot output directory not given (--plot-output)')
        if not options.json_output:
            parser.error('Json output file not given (--json)')
        main(options)
    except optparse.OptionError as exc:
        raise RuntimeError(f"Invalid arguments: {args}") from exc
