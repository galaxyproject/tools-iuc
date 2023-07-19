import logging
import os
from typing import List

from data import Dataset

from hydra.core.config_store import ConfigStore

import matplotlib
import matplotlib.pyplot as plt

from omegaconf import DictConfig

import pandas as pd

import seaborn as sns


logger = logging.getLogger(__name__)

cs = ConfigStore.instance()


def melt_data_metadata_2df(compartment_df: pd.DataFrame,
                           metadata_co_df: pd.DataFrame) -> pd.DataFrame:
    """
    - merges quantitative data and metadata
    - melts the merged df to obtain df:
         "timenum", "condition", "metabolite" "Fractional Contribution (%)"
    - quantitative values are multiplied by 100, as expressed in %
    """
    compartment_df = compartment_df.T
    compartment_df['name_to_plot'] = compartment_df.index
    compartment_df = pd.merge(compartment_df, metadata_co_df,
                              on='name_to_plot')
    compartment_df = compartment_df.drop(columns=['name_to_plot',
                                                  'timepoint',
                                                  'short_comp',
                                                  'original_name'])
    melted_df = pd.melt(compartment_df,
                        id_vars=['timenum', 'condition'],
                        var_name="metabolite",
                        value_name="Fractional Contribution (%)")
    melted_df["Fractional Contribution (%)"] = \
        melted_df["Fractional Contribution (%)"] * 100
    return melted_df


def nested_dict__2_list(nested_dict) -> List[str]:
    result = []
    for metabolites_list in nested_dict.values():
        for metabolite in metabolites_list:
            result.append(metabolite)
    return result


def metabolite_df__mean_and_sd(
        one_metabolite_df: pd.DataFrame) -> pd.DataFrame:
    """
    input: dataframe by one metabolite:
      "timenum", "condition", "metabolite" "Fractional Contribution (%)"
    returns :
        dataframe by one metabolite
                condition  timenum    mean    sd    metabolite
    108        Control     0  0.000893  0.002611
    111        Control     1  0.236453  0.023246
    ...
    123        Control    24  0.854101  0.055241
    126  L-Cycloserine     0  0.010083  0.003465
    129  L-Cycloserine     1  0.259570  0.008602
    ...
    141  L-Cycloserine    24  0.815613  0.050756
    """
    assert len(one_metabolite_df["metabolite"].unique()) == 1
    df = one_metabolite_df.copy()
    df = df.drop_duplicates()
    # by default both std and mean in pandas ignore NaN
    mean_df = df.groupby(["condition", "timenum", "metabolite"])[
        "Fractional Contribution (%)"].mean().reset_index(name="mean")
    # std by pandas : ddof=0 to have same result as with numpy std
    std_df = df.groupby(["condition", "timenum", "metabolite"])[
        "Fractional Contribution (%)"].std(ddof=0).reset_index(name="sd")

    one_metabolite_result = mean_df.merge(std_df, how='inner',
                                          on=["condition", "timenum",
                                              "metabolite"])
    return one_metabolite_result


def add_mean_and_sd__df(
        metabolites_selected_df: pd.DataFrame) -> pd.DataFrame:
    tmp_list: List[pd.DataFrame] = list()
    for metabolite_i in set(metabolites_selected_df["metabolite"]):
        one_metabolite_df = metabolites_selected_df[
            metabolites_selected_df["metabolite"] == metabolite_i]
        tmp_df = metabolite_df__mean_and_sd(one_metabolite_df)
        tmp_list.append(tmp_df)
    result_df = pd.concat(tmp_list, axis=0)
    return result_df


def plot_one_metabolite(metabolite_name: str,
                        metabolite_df: pd.DataFrame,
                        mean_and_sd_metabolite_df: pd.DataFrame,
                        color_lines_by: str,
                        palette_option: dict,
                        xaxis_title: str,
                        time_ticks: List[float],
                        alpha_conf: float,
                        axs_z: matplotlib.axes
                        ) -> matplotlib.axes:
    """
    Often single metabolite is plotted, as configured by default.

    If several metabolites in input dataframes,
    they will be plotted sharing axes in same plane.
    """
    sns.lineplot(
        ax=axs_z,  # starts in 1 (after the empty row index 0)
        x="timenum",
        y="Fractional Contribution (%)",
        hue=color_lines_by,
        style="condition",
        err_style=None,
        alpha=alpha_conf,
        linewidth=4.5,
        palette=palette_option[color_lines_by],
        zorder=1,
        data=metabolite_df,

        legend=True,
    )
    axs_z.set_xticks([float(i) for i in time_ticks])

    axs_z.scatter(
        mean_and_sd_metabolite_df["timenum"],
        mean_and_sd_metabolite_df["mean"], s=23,
        facecolors="none", edgecolors="black"
    )
    axs_z.errorbar(
        mean_and_sd_metabolite_df["timenum"],
        mean_and_sd_metabolite_df["mean"],
        yerr=mean_and_sd_metabolite_df["sd"],
        fmt="none",
        capsize=3.5,
        ecolor="black",
        zorder=2
    )
    axs_z.set_ylabel("Fractional Contribution (%)", size=19),
    axs_z.set(xlabel=xaxis_title)
    axs_z.set(title=metabolite_name)
    axs_z.legend(loc="upper left",
                 bbox_to_anchor=(1, 1), ncol=1,
                 frameon=False)
    return axs_z


def save_line_plot_no_grid(melted_compartment_df: pd.DataFrame,
                           metabolites_compartment_dict: dict,
                           width_subplot: float,
                           height_subplot: float,
                           color_lines_by: str,
                           palette_option: dict,
                           xaxis_title: str,
                           alpha_conf: float,
                           out_file_elements: List[str],
                           out_plot_dir: str) -> None:
    """Save each metabolite in independent plots"""
    time_ticks = melted_compartment_df['timenum'].unique()
    metabolites_selected_df = \
        melted_compartment_df.loc[melted_compartment_df["metabolite"].isin(
            nested_dict__2_list(metabolites_compartment_dict)), :].copy()

    mean_and_sd_df = add_mean_and_sd__df(metabolites_selected_df)

    sns.set_style(
        {"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"}
    )
    plt.rcParams.update({"font.size": 22})

    for z in range(len(metabolites_compartment_dict)):
        fig, axs_unique = plt.subplots(
            1, 1,
            sharey=False,
            figsize=(width_subplot, height_subplot)
        )

        metabolite_df = metabolites_selected_df.loc[
            metabolites_selected_df["metabolite"].isin(
                metabolites_compartment_dict[z])]

        mean_and_sd_metabolite_df = mean_and_sd_df.loc[
            mean_and_sd_df["metabolite"].isin(
                metabolites_compartment_dict[z])]
        metabolite_name = str(metabolites_compartment_dict[z][0])
        plot_one_metabolite(  # implicitly returns the axs object
            metabolite_name,
            metabolite_df,
            mean_and_sd_metabolite_df,
            color_lines_by,
            palette_option,
            xaxis_title,
            time_ticks,
            alpha_conf,
            axs_unique
        )
        out_file_elements_plus_met = out_file_elements + [metabolite_name]
        output_path = os.path.join(
            out_plot_dir, f'{"-".join(out_file_elements_plus_met)}.pdf')
        fig.savefig(output_path, format="pdf", bbox_inches='tight')
        plt.close()
    # end for
    logger.info(f"Saved mean enrichment line plots in {out_plot_dir}")


def save_line_plot_as_grid(melted_compartment_df: pd.DataFrame,
                           metabolites_compartment_dict: dict,
                           width_subplot: float,
                           height_subplot: float,
                           color_lines_by: str,
                           palette_option: dict,
                           xaxis_title: str,
                           alpha_conf: float,
                           out_file_elements: List[str],
                           out_plot_dir) -> None:
    """
    constructs and saves the grid plot to pdf
    """
    time_ticks = melted_compartment_df['timenum'].unique()
    metabolites_selected_df = \
        melted_compartment_df.loc[melted_compartment_df["metabolite"].isin(
            nested_dict__2_list(metabolites_compartment_dict)), :].copy()

    mean_and_sd_df = add_mean_and_sd__df(metabolites_selected_df)

    corrector_factor = 0.7
    total_height_grid = height_subplot * len(metabolites_compartment_dict) + \
        (corrector_factor * 7 * len(metabolites_compartment_dict))

    sns.set_style(
        {"font.family": "sans-serif", "font.sans-serif": "Liberation Sans"}
    )
    plt.rcParams.update({"font.size": 22})
    fig, axs = plt.subplots(
        len(metabolites_compartment_dict),
        1,
        sharey=False,
        figsize=(width_subplot, total_height_grid)
    )
    for z in range(len(metabolites_compartment_dict)):
        metabolite_df = metabolites_selected_df.loc[
            metabolites_selected_df["metabolite"].isin(
                metabolites_compartment_dict[z])]

        mean_and_sd_metabolite_df = mean_and_sd_df.loc[
            mean_and_sd_df["metabolite"].isin(
                metabolites_compartment_dict[z])]

        metabolite_name = str(metabolites_compartment_dict[z][0])
        axs[z] = plot_one_metabolite(
            metabolite_name,
            metabolite_df,
            mean_and_sd_metabolite_df,
            color_lines_by,
            palette_option,
            xaxis_title,
            time_ticks,
            alpha_conf,
            axs[z]
        )
    plt.subplots_adjust(hspace=corrector_factor)

    output_path = os.path.join(out_plot_dir,
                               f'{"-".join(out_file_elements)}.pdf')
    fig.savefig(output_path, format="pdf", bbox_inches='tight')
    plt.close
    logger.info(f"Saved mean enrichment line plots in {output_path}")


def give_colors_by_metabolite(cfg: DictConfig,
                              metabolites_numbered_dict) -> dict:
    handycolors = ["rosybrown", "lightcoral", "brown", "firebrick",
                   "tomato", "coral", "sienna", "darkorange", "peru",
                   "darkgoldenrod", "gold", "darkkhaki", "olive",
                   "yellowgreen", "limegreen", "green", "lightseagreen",
                   "mediumturquoise", "darkcyan", "teal", "cadetblue",
                   "slategrey", "steelblue", "navy", "darkslateblue",
                   "blueviolet",
                   "darkochid", "purple", "mediumvioletred", "crimson"]

    colors_dict = dict()

    if cfg.analysis.method.palette_metabolite == "auto_multi_color":
        tmp = set()
        for co in metabolites_numbered_dict.keys():
            for k in metabolites_numbered_dict[co].keys():
                tmp.update(set(metabolites_numbered_dict[co][k]))
        metabolites = sorted(list(tmp))
        if len(metabolites) <= 12:
            # default Paired palette for coloring individual metabolites
            palettecols = sns.color_palette("Paired", 12)
            for i in range(len(metabolites)):
                colors_dict[metabolites[i]] = palettecols[i]
        else:  # more than 12 colors obliges to set them manually
            for i in range(len(metabolites)):
                colors_dict[metabolites[i]] = handycolors[i]
    else:  # argument_color_metabolites is a csv file
        try:
            file_containing_colors = cfg.analysis.method.palette_metabolite
            df = pd.read_csv(file_containing_colors, header=0)
            for i, row in df.iterrows():
                metabolite = df.iloc[i, 0]  # first column metabolite
                color = df.iloc[i, 1]  # second column is color
                colors_dict[metabolite] = color
        except Exception as e:
            logger.info(e, f"\n could not assign color, wrong csv file: \
                   {cfg.analysis.method.palette_metabolite}")
            colors_dict = None

    return colors_dict


def give_colors_by_option(cfg: DictConfig,
                          metabolites_numbered_dict) -> dict:
    """
    if option color_lines_by = metabolite, returns a dictionary of colors;
    otherwise color_lines_by = condition, returns a string (palette name)
    """
    assert cfg.analysis.method.color_lines_by in ["metabolite", "condition"]
    if cfg.analysis.method.color_lines_by == "metabolite":
        colors_dict: dict = give_colors_by_metabolite(
            cfg, metabolites_numbered_dict)
        palette_option = {
            "metabolite": colors_dict
        }
    else:
        try:
            palette_condition: str = cfg.analysis.method.palette_condition
        except ValueError:
            palette_condition: str = "paired"
        palette_option = {
            "condition": palette_condition
        }
    return palette_option


def line_plot_by_compartment(dataset: Dataset,
                             conditions_leveled: List[str],
                             out_plot_dir: str,
                             metabolites_numbered_dict,
                             cfg: DictConfig) -> None:
    """ calls function to construct and save plot """
    metadata_df = dataset.metadata_df
    compartments = list(metadata_df['short_comp'].unique())
    width_subplot = cfg.analysis.width_subplot
    height_subplot = cfg.analysis.method.height_subplot
    xaxis_title = cfg.analysis.method.xaxis_title
    color_lines_by = cfg.analysis.method.color_lines_by
    palette_option = give_colors_by_option(cfg, metabolites_numbered_dict)
    alpha_conf = cfg.analysis.method.alpha

    for co in compartments:
        metadata_co_df = metadata_df.loc[metadata_df['short_comp'] == co, :]
        compartment_df = dataset.compartmentalized_dfs["mean_enrichment"][co]

        melted_co_df = melt_data_metadata_2df(compartment_df, metadata_co_df)
        melted_co_df["condition"] = pd.Categorical(
            melted_co_df["condition"], conditions_leveled)
        metabolites_compartment_dict = metabolites_numbered_dict[co]

        out_file_elements = ["mean_enrichment", co]

        save_line_plot_no_grid(melted_co_df, metabolites_compartment_dict,
                               width_subplot,
                               height_subplot,
                               color_lines_by,
                               palette_option,
                               xaxis_title,
                               alpha_conf,
                               out_file_elements, out_plot_dir)

        if (cfg.analysis.method.as_grid is not None) and \
                cfg.analysis.method.as_grid:
            save_line_plot_as_grid(
                melted_co_df,
                metabolites_compartment_dict,
                width_subplot,
                height_subplot,
                color_lines_by,
                palette_option,
                xaxis_title,
                alpha_conf,
                out_file_elements, out_plot_dir)


def generate_metabolites_numbered_dict(cfg: DictConfig,
                                       metabolites: List[str]) -> dict:
    """
    Helps the user in the case of needing two metabolites in a same plane,
    by using 'plot_grouped_by_dict' option, example:
    { 0: ['Pyr', 'Cit'], 1: ['Asn', 'Asp']}
    will result in one plot by couple of metabolites, totalling 2 plots.
    If option 'plot_grouped_by_dict' not specified, one single metabolite
    per numeric key is set.
    """
    if cfg.analysis.method.plot_grouped_by_dict is not None:
        metabolites_numbered_dict = cfg.analysis.method.plot_grouped_by_dict

    else:
        try:
            tmp = dict()
            for compartment in metabolites.keys():
                tmp[compartment] = dict()
                for i, m in enumerate(metabolites[compartment]):
                    tmp[compartment][i] = [m]
            metabolites_numbered_dict = tmp
        except KeyError:
            logger.info("run_mean_enrichment_line_plot: \
               No metabolites for plotting in your config file")

    return metabolites_numbered_dict


def run_mean_enrichment_line_plot(dataset: Dataset,
                                  out_plot_dir: str,
                                  cfg: DictConfig) -> None:
    metabolites: dict = (
        cfg.analysis.metabolites
    )  # will define which metabolites are plotted
    conditions_leveled = cfg.analysis.dataset.conditions

    metabolites_numbered_dict = generate_metabolites_numbered_dict(
        cfg, metabolites
    )

    line_plot_by_compartment(dataset,
                             conditions_leveled,
                             out_plot_dir,
                             metabolites_numbered_dict,
                             cfg)
