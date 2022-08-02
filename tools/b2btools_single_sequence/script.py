import optparse
import os.path
import unicodedata
import re
import numpy as np
import pandas as pd
from b2bTools import SingleSeq
import matplotlib.pyplot as plt


def slugify(value):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode'. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


def check_min_max(predicted_values, former_min, former_max):
    seq_max = max(predicted_values)
    seq_min = min(predicted_values)
    if seq_max + 0.1 > former_max and not np.isnan(seq_max) and not np.isinf(seq_max):
        former_max = seq_max + 0.1
    if seq_min - 0.1 < former_min and not np.isnan(seq_min) and not np.isinf(seq_min):
        former_min = seq_min - 0.1
    return former_min, former_max


def plot_prediction(prediction_name, highlighting_regions, predicted_values, seq_name):
    thresholds_dict = {'backbone': {'membrane spanning': [1., 1.5],
                                    'rigid': [0.8, 1.],
                                    'context-dependent': [0.69, 0.8],
                                    'flexible': [-1.0, 0.69]},
                       'earlyFolding': {'early folds': [0.169, 2.], 'late folds': [-1., 0.169]},
                       'disoMine': {'ordered': [-1., 0.5], 'disordered': [0.5, 2.]},
                       }
    ordered_regions_dict = {'backbone': ['flexible', 'context-dependent', 'rigid', 'membrane spanning'],
                            'earlyFolding': ['late folds', 'early folds'],
                            'disoMine': ['ordered', 'disordered'],
                            }
    colors = ['yellow', 'orange', 'pink', 'red']
    ranges_dict = {
        'backbone': [-0.2, 1.2],
        'sidechain': [-0.2, 1.2],
        'ppII': [-0.2, 1.2],
        'earlyFolding': [-0.2, 1.2],
        'disoMine': [-0.2, 1.2],
        'agmata': [-0.2, 1.2],
        'helix': [-1., 1.],
        'sheet': [-1., 1.],
        'coil': [-1., 1.],
    }
    fig, ax = plt.subplots(1, 1)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    ax.set_title(prediction_name + ' ' + 'prediction')
    min_value, max_value = ranges_dict[prediction_name]
    if seq_name == 'all':
        max_len = 0
        for seq in predicted_values.keys():
            predictions = predicted_values[seq]
            min_value, max_value = check_min_max(predictions, min_value, max_value)
            ax.plot(range(len(predictions)), predictions, label=seq)
            if len(predictions) > max_len:
                max_len = len(predictions)
            ax.set_xlim([0, max_len - 1])
    else:
        predictions = predicted_values
        min_value, max_value = check_min_max(predictions, min_value, max_value)
        ax.plot(range(len(predictions)), predictions, label=seq_name)
        ax.set_xlim([0, len(predictions) - 1])
    legend_lines = plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fancybox=True, shadow=True)
    ax.add_artist(legend_lines)
    # Define regions
    if highlighting_regions:
        if prediction_name in ordered_regions_dict.keys():
            for i, prediction in enumerate(ordered_regions_dict[prediction_name]):
                lower = thresholds_dict[prediction_name][prediction][0]
                upper = thresholds_dict[prediction_name][prediction][1]
                color = colors[i]
                ax.axhspan(lower, upper, alpha=0.3, color=color, label=prediction)
            included_in_regions_legend = list(reversed(
                [prediction for prediction in ordered_regions_dict[prediction_name]]))  # to sort it "from up to low"
            # Get handles and labels
            handles, labels = plt.gca().get_legend_handles_labels()
            handles_dict = {label: handles[idx] for idx, label in enumerate(labels)}
            # Add legend for regions, if available
            region_legend = ax.legend([handles_dict[region] for region in included_in_regions_legend],
                                      [region for region in included_in_regions_legend], fancybox=True, shadow=True,
                                      loc='lower left', bbox_to_anchor=(1.04, 0))
            ax.add_artist(region_legend)
    ax.set_ylim([min_value, max_value])
    ax.set_xlabel('residue index')
    ax.set_ylabel('prediction values')
    ax.grid(axis='y')
    plt.savefig(os.path.join(options.plot_output, "{0}_{1}.png".format(slugify(seq_name), prediction_name)), bbox_inches="tight")
    plt.close()


def df_dict_to_dict_of_values(df_dict, predictor):
    results_dict = {}
    for seq in df_dict.keys():
        df = pd.read_csv(df_dict[seq], sep='\t')
        results_dict[seq] = df[predictor]
    return results_dict


def main(options):
    single_seq = SingleSeq(options.input_fasta)
    b2b_tools = []
    if options.dynamine:
        b2b_tools.append('dynamine')
    if options.disomine:
        b2b_tools.append('disomine')
    if options.efoldmine:
        b2b_tools.append('efoldmine')
    if options.agmata:
        b2b_tools.append('agmata')

    single_seq.predict(b2b_tools)
    predictions = single_seq.get_all_predictions()
    results_json = single_seq.get_all_predictions_json('all')
    with open(options.json_output, 'w') as f:
        f.write(results_json)
    first_sequence_key = next(iter(predictions))
    prediction_keys = predictions[first_sequence_key].keys()
    df_dictionary = {}
    for sequence_key, sequence_predictions in predictions.items():
        residues = sequence_predictions['seq']
        residues_count = len(residues)
        sequence_df = pd.DataFrame(columns=prediction_keys, index=range(residues_count))
        sequence_df.index.name = 'residue_index'
        for predictor in prediction_keys:
            sequence_df[predictor] = sequence_predictions[predictor]
        sequence_df = sequence_df.rename(columns={"seq": "residue"})
        sequence_df = sequence_df.round(decimals=2)
        filename = f'{options.output}/{slugify(sequence_key)}.tsv'
        df_dictionary[sequence_key] = filename
        sequence_df.to_csv(filename, sep="\t")
        # Plot each individual plot (compatible with plot all)
        if options.plot:
            for predictor in prediction_keys:
                if predictor != 'seq':
                    plot_prediction(prediction_name=predictor, highlighting_regions=True,
                                    predicted_values=sequence_predictions[predictor], seq_name=sequence_key)
    # Plot all together (compatible with plot individual)
    if options.plot_all:
        for predictor in prediction_keys:
            if predictor != 'seq':
                results_dictionary = df_dict_to_dict_of_values(df_dict=df_dictionary, predictor=predictor)
                plot_prediction(prediction_name=predictor, highlighting_regions=True,
                                predicted_values=results_dictionary, seq_name='all')


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("--dynamine", action="store_true", default=False)
    parser.add_option("--disomine", action="store_true", default=False)
    parser.add_option("--efoldmine", action="store_true", default=False)
    parser.add_option("--agmata", action="store_true", default=False)
    parser.add_option("--file", dest="input_fasta", default=False)
    parser.add_option("--output", dest="output", default=False)
    parser.add_option("--plot-output", dest="plot_output", default=False)

    parser.add_option("--json", dest="json_output", default=False)
    parser.add_option("--plot", action="store_true", default=False)
    parser.add_option("--plot_all", action="store_true", default=False)
    parser.add_option("--highlight", action="store_true", default=False)
    options, _args = parser.parse_args()
    main(options)
