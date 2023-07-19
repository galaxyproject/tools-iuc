#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
import logging
import os
from typing import Dict, List

from data import Dataset

from hydra.core.config_store import ConfigStore

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from omegaconf import DictConfig

import pandas as pd

import seaborn as sns


logger = logging.getLogger(__name__)

cs = ConfigStore.instance()


def isotopologue_proportions_2piled_df(
        compartment_df: pd.DataFrame,
        metada_df: pd.DataFrame) -> pd.DataFrame:
    """
    melt the compartment_df, several steps required:
    - transpose compratment_df
    - combine compratment_df and metadata
    - melt
    - multiply isotopologues by 100 (are given as proportions, must be %
        for plotting)
    - set
    example output:
        timenum   condition  isotopologue_name  Isotopologue Contribution (%)
            0    control     L-Phenylalanine     0.01
    """
    combined_isos_metadata_df = compartment_df.T

    combined_isos_metadata_df[
        'name_to_plot'] = combined_isos_metadata_df.index
    combined_isos_metadata_df = pd.merge(
        combined_isos_metadata_df, metada_df, on='name_to_plot')

    combined_isos_metadata_df = combined_isos_metadata_df.drop(
        columns=['short_comp', 'original_name', 'name_to_plot', 'timepoint'])
    piled_df = pd.melt(combined_isos_metadata_df,
                       id_vars=['timenum', 'condition'],
                       var_name="isotopologue_name",
                       value_name="Isotopologue Contribution (%)")

    piled_df['timenum'] = piled_df['timenum'].astype(str)
    piled_df['Isotopologue Contribution (%)'] = \
        piled_df['Isotopologue Contribution (%)'] * 100
    return piled_df


def massage_isotopologues(piled_df) -> pd.DataFrame:
    """
    returns dataframe splitting metabolite and m+x into two separate columns
    and also correcting weird values
    """
    tmp_df = piled_df['isotopologue_name'].str.split("_m+",
                                                     expand=True, regex=False)
    tmp_df.rename(columns={0: 'name', 1: 'm+x'}, inplace=True)
    piled_df["metabolite"] = tmp_df["name"]
    piled_df["m+x"] = tmp_df["m+x"]
    piled_df["m+x"] = "m+" + tmp_df["m+x"].astype(str)

    # dealing with weird values: bigger than 100 and less than 0 :
    piled_df.loc[
        piled_df["Isotopologue Contribution (%)"] > 100,
        "Isotopologue Contribution (%)"
    ] = 100

    piled_df.loc[
        piled_df["Isotopologue Contribution (%)"] < 0,
        "Isotopologue Contribution (%)"
    ] = 0

    return piled_df


def prepare_means_replicates(piled_df, metaboli_selected) -> Dict:
    """
    returns a dictionary of dataframes, keys are metabolites
    """
    dfcopy = piled_df.copy()
    # instead groupby isotopologue_name, using m+x and metabolite works better
    dfcopy = dfcopy.groupby(
        ["condition", "metabolite", "m+x", "timenum"]) \
        .mean("Isotopologue Contribution %")  # df.mean skips nan by default
    dfcopy = dfcopy.reset_index()

    dfs_dict = dict()
    for i in metaboli_selected:
        tmp = dfcopy.loc[dfcopy["metabolite"] == i, ].reset_index(drop=True)
        # set m+x as numeric to avoid any bad reordering of stacked m+x
        tmp["m+x"] = tmp["m+x"].str.split("m+", regex=False).str[1]
        tmp["m+x"] = tmp["m+x"].astype(int)

        dfs_dict[i] = tmp
    return dfs_dict


def add_combined_conditime(dfs_dict: Dict,
                           combined_tc_levels: List[str]) -> Dict:
    """
    add column 'time_and_condition' to each  metabolite dataframe in Dictio
    """
    for metab in dfs_dict.keys():
        dfs_dict[metab]["time_and_condition"] = \
            dfs_dict[metab]["timenum"] + " : " + dfs_dict[metab]["condition"]

        dfs_dict[metab]["time_and_condition"] = pd.Categorical(
            dfs_dict[metab]["time_and_condition"],
            combined_tc_levels)
    return dfs_dict


def add_categorical_time(dfs_dict, levelstime_str) -> Dict:
    """
    existing column 'timenum' as categorical, with specific order.
    Each metabolite dataframe in Dictio is processed.
    """
    for metab in dfs_dict.keys():
        dfs_dict[metab]["timenum"] = pd.Categorical(
            dfs_dict[metab]["timenum"],
            levelstime_str)
    return dfs_dict


def give_colors_carbons(nb_of_carbons: int) -> Dict:
    """
    currently 30 carbons (30 colors) supported
    """
    color_d = dict()
    # color_d[0] = "lightgray"  # m+0
    color_d[0] = "#410257"
    # set colors m+1 to m+8 from Spectral palette,
    # with custom spaced selected colors (validated)
    spectralPal = sns.color_palette("Spectral", 30)
    color_d[1] = spectralPal[29]
    color_d[2] = spectralPal[26]
    color_d[3] = spectralPal[21]
    color_d[4] = spectralPal[17]
    color_d[5] = spectralPal[10]
    color_d[6] = spectralPal[6]
    color_d[7] = spectralPal[3]
    color_d[8] = spectralPal[0]
    # rest of the colors from tab20b palette
    added_pal = sns.color_palette("tab20b", 20)
    i = 9
    j = 19
    while i <= nb_of_carbons:
        color_d[i] = added_pal[j]
        j = j - 1
        i += 1

    return color_d


def plot_one_metabolite(
        metabolite_name: str,
        one_metabolite_df: pd.DataFrame,
        x_to_plot: str,
        numbers_size: float,
        x_ticks_text_size: float,
        x_ticks_text_tilt: float,
        colors_isotopologues_dict: dict,
        curr_ax: matplotlib.axes
):
    curr_ax.set_title(metabolite_name)
    sns.histplot(
        ax=curr_ax,
        data=one_metabolite_df,
        x=x_to_plot,
        # Use the value variable here to turn histogram counts into
        # weighted values.
        weights="Isotopologue Contribution (%)",
        hue="m+x",
        multiple="stack",
        palette=colors_isotopologues_dict,
        # Add  borders to the bars.
        edgecolor="black",
        # Shrink the bars a bit, so they don't touch.
        shrink=0.85,
        alpha=1,
        legend=False,
    )

    for xtick in curr_ax.get_xticklabels():
        if xtick.get_text().endswith("xemptyspace"):
            xtick.set_color("white")
        else:
            xtick.set_color("black")

    curr_ax.tick_params(axis="x",
                        labelrotation=x_ticks_text_tilt,
                        labelsize=x_ticks_text_size)

    curr_ax.tick_params(axis="y", length=3,
                        labelsize=19)
    curr_ax.set_ylim([0, 100])

    # #  * * Inner text * *
    if numbers_size <= 0:
        numbers_alpha = 0  # to make it invisible
    else:
        numbers_alpha = 1
    # for defining the color of inner text in bar when M0
    rgba_eq_hex_410257 = \
        (0.2549019607843137, 0.00784313725490196, 0.3411764705882353, 1.0)

    for bar in curr_ax.patches:
        inner_text_color = "black"
        here_rgba = bar.get_facecolor()
        if here_rgba == rgba_eq_hex_410257:
            inner_text_color = "white"
        thebarvalue = round(bar.get_height(), 1)
        if thebarvalue >= 100:
            thebarvalue = 100  # no decimals if 100
        if round(bar.get_height(), 1) >= 4:
            curr_ax.text(
                # Put the text in the middle of each bar. get_x returns t
                # he start, so we add half the width to get to the middle.
                bar.get_x() + bar.get_width() / 2,
                # Vertically, add the height of the bar to the start of
                # the bar, along with the offset.
                (bar.get_height() / 2) + (bar.get_y()) + 2,  #
                # This is actual value we'll show.
                thebarvalue,
                # Center the labels and style them a bit.
                ha="center",
                color=inner_text_color,
                size=numbers_size,
                alpha=numbers_alpha
            )  # end curr_ax.text
        else:
            continue
        # end if round
    # end for bar; end Inner text

    # #  * * y labels ticks and labels organize * *
    curr_ax.set_ylabel("Isotopologue Contribution (%)", size=19)
    curr_ax.xaxis.set_tick_params(length=0)  # no need of x ticks
    curr_ax.set_xlabel("", size=13)
    curr_ax.invert_yaxis()  # invert y, step 1    #

    ylabels = list(curr_ax.get_yticks())
    curr_ax.yaxis.set_major_locator(mticker.FixedLocator(ylabels))
    curr_ax.set_yticklabels(
        [100 - int(i) for i in ylabels])  # invert y, step2

    return curr_ax


def complex_stacked_plot_as_grid(
        metaboli_selected: List[str], dfs_dict: Dict,
        out_plot_dir: str, output_file_elements: List[str],
        cfg: DictConfig,
        xlab_text: bool,
        x_to_plot: str, x_ticks_text_tilt: int) -> None:
    """
    Using the isotopologues proportions, generates stacked barplots
    per metabolite, all arranged in a single pdf file.
    A legend is produced separately also as pdf file.
    """
    width_each_stack = cfg.analysis.width_each_stack
    numbers_size = cfg.analysis.inner_numbers_size
    height_each_stack = cfg.analysis.method.height_each_stack

    x_ticks_text_size = cfg.analysis.method.x_ticks_text_size
    colors_isotopologues_dict = give_colors_carbons(
        cfg.analysis.method.max_nb_carbons_possible)

    corrector_factor = 0.7
    total_height_grid = height_each_stack * len(metaboli_selected) + \
        (corrector_factor * 7 * len(metaboli_selected))

    f, axs = plt.subplots(len(metaboli_selected), 1,
                          sharey=cfg.analysis.method.sharey,
                          figsize=(width_each_stack,
                                   total_height_grid))
    plt.rcParams.update({"font.size": 20})

    for z in range(len(metaboli_selected)):
        one_metabolite_name: str = metaboli_selected[z]
        one_metabolite_df = dfs_dict[one_metabolite_name]
        try:
            axs[z] = plot_one_metabolite(
                one_metabolite_name,
                one_metabolite_df,
                x_to_plot,
                numbers_size,
                x_ticks_text_size,
                x_ticks_text_tilt,
                colors_isotopologues_dict, axs[z]
            )
        except ValueError:
            if one_metabolite_df[
                    "Isotopologue Contribution (%)"].isnull().all():
                axs[z].set_title(one_metabolite_name + "\n" +
                                 "NaN values only")
            pass
        except Exception as e:
            print(e)
            pass
    # end for z
    # when panel without x labels
    if not xlab_text:
        for ax in axs:
            ax.get_xaxis().set_ticks([])

    f.subplots_adjust(hspace=corrector_factor, top=0.85,
                      bottom=0.26, left=0.15, right=0.99)

    file_elems_str = "-".join(output_file_elements)
    output_path = os.path.join(
        out_plot_dir, f"Isotopologues_{file_elems_str}.pdf"
    )
    f.savefig(output_path,
              bbox_inches="tight", format="pdf")
    plt.close()
    logger.info(f"Saved isotopologue stacked barplots in {output_path}")


def stacked_plot_no_grid(
        metaboli_selected: List[str], dfs_dict: Dict,
        out_plot_dir: str, output_file_elements: List[str],
        cfg: DictConfig,
        xlab_text: bool,
        x_to_plot: str, x_ticks_text_tilt: int) -> None:
    """
    Using the isotopologues proportions, generates stacked barplots
    per metabolite, each one in independent pdf file
    """
    width_each_stack = cfg.analysis.width_each_stack
    numbers_size = cfg.analysis.inner_numbers_size
    height_each_stack = cfg.analysis.method.height_each_stack
    x_ticks_text_size = cfg.analysis.method.x_ticks_text_size
    colors_isotopologues_dict = give_colors_carbons(
        cfg.analysis.method.max_nb_carbons_possible)

    plt.rcParams.update({"font.size": 20})

    for z in range(len(metaboli_selected)):
        f, axs_z = plt.subplots(1, 1,
                                sharey=cfg.analysis.method.sharey,
                                figsize=(width_each_stack,
                                         height_each_stack))
        one_metabolite_name: str = metaboli_selected[z]
        one_metabolite_df = dfs_dict[one_metabolite_name]
        try:
            axs_z = plot_one_metabolite(
                one_metabolite_name,
                one_metabolite_df,
                x_to_plot,
                numbers_size,
                x_ticks_text_size,
                x_ticks_text_tilt,
                colors_isotopologues_dict, axs_z
            )
        except ValueError:
            if one_metabolite_df[
                    "Isotopologue Contribution (%)"].isnull().all():
                axs_z.set_title(one_metabolite_name + "\n" +
                                "NaN values only")
            pass
        except Exception as e:
            print(e)
            pass
        # when panel without x labels
        if not xlab_text:
            axs_z.get_xaxis().set_ticks([])

        file_elems_str = "-".join(
            output_file_elements + [one_metabolite_name])
        output_path = os.path.join(
            out_plot_dir, f"Isotopologues_{file_elems_str}.pdf"
        )
        f.savefig(output_path,
                  bbox_inches="tight", format="pdf")
        plt.close()
    # end for
    logger.info(f"Saved isotopologue stacked barplots in {output_path}")


def add_xemptyspace_tolabs(conditions, time_levels_list):
    """
    adds an 'empty space' between each timepoint in the metabolite plot,
    conditions are kept next to each other in comparative aspect.
    Note : If willing to _not_ having conditions in comparative aspect
      but each condition in separate pdf instead,
      see .split_plots_by_condition attribute.
    """
    conditions.extend(["xemptyspace"])  # to add space among time categories
    combined_tc_levels = list()
    tmp = conditions.copy()
    conditions = list()  # warranty uniqueness
    for i in tmp:
        if i not in conditions:
            conditions.append(i)
    for x in time_levels_list:
        for y in conditions:
            if y == "xemptyspace":
                combined_tc_levels.append(str(x) + "xemptyspace")
            else:
                combined_tc_levels.append(f'{x} : {y}')
    return conditions, combined_tc_levels


def time_plus_condi_labs(conditions, time_levels_list):
    combined_tc_levels = list()
    for x in time_levels_list:
        for y in conditions:
            combined_tc_levels.append(f'{x} : {y}')
    return combined_tc_levels


def give_labels_to_colors_dict(colors_isotopologues_dict) -> Dict:
    tmp = dict()
    for k in colors_isotopologues_dict.keys():
        tmp["m+" + str(k)] = colors_isotopologues_dict[k]
    return tmp


def isotopologues_types_create_legend(cfg) -> matplotlib.figure.Figure:
    plt.figure(figsize=(4, cfg.analysis.method.max_nb_carbons_possible * 0.6))
    colors_isotopologues_dict = give_colors_carbons(
        cfg.analysis.method.max_nb_carbons_possible)
    colors_isotopologues_dict_labeled = give_labels_to_colors_dict(
        colors_isotopologues_dict)
    myhandless = []
    for c in colors_isotopologues_dict_labeled.keys():
        paobj = mpatches.Patch(facecolor=colors_isotopologues_dict_labeled[c],
                               label=c, edgecolor="black")
        myhandless.append(paobj)
    plt.legend(handles=myhandless, labelspacing=0.01)
    plt.axis("off")
    return plt


def run_isotopologue_proportions_plot(dataset: Dataset,
                                      out_plot_dir, cfg: DictConfig) -> None:
    metadata_df = dataset.metadata_df
    timepoints = cfg.analysis.timepoints
    metabolites = cfg.analysis.metabolites
    conditions = cfg.analysis.dataset.conditions
    x_ticks_text_tilt_fixed: dict = {  # usr def tilt bad result, let fixed
        'time_only': 0, 'time_and_condition': 90
    }
    time_levels_list: List[str] = [
        str(i) for i in sorted(metadata_df['timenum'].unique())]

    compartments = list(metadata_df['short_comp'].unique())

    for compartment in compartments:
        metadata_compartment_df: pd.DataFrame = \
            metadata_df.loc[metadata_df["short_comp"] == compartment, :]
        compartment_df = dataset.compartmentalized_dfs[
            "isotopologue_proportions"][compartment]

        # metadata, isotopologues and time of interest
        time_metadata_df = metadata_compartment_df.loc[
                           metadata_compartment_df["timepoint"].isin(
                               timepoints), :]
        time_compartment_df = compartment_df[time_metadata_df["name_to_plot"]]
        # note that pandas automatically transforms any 99.9% in decimal 0.999
        piled_df = isotopologue_proportions_2piled_df(time_compartment_df,
                                                      time_metadata_df)
        # values now are in %
        piled_df = massage_isotopologues(piled_df)
        compartment_metabolites = metabolites[compartment]
        dfs_dict = prepare_means_replicates(piled_df, compartment_metabolites)

        if cfg.analysis.method.split_plots_by_condition or \
                (len(conditions) == 1):
            plt.set_loglevel('WARNING')  # disable INFO matplotlib
            for condition in conditions:
                cond_time_metadata = time_metadata_df.loc[
                                     time_metadata_df[
                                         'condition'] == condition, :]
                condition_df = time_compartment_df[
                    cond_time_metadata['name_to_plot']]
                piled_df = isotopologue_proportions_2piled_df(
                    condition_df, cond_time_metadata)
                piled_df = massage_isotopologues(piled_df)
                dfs_dict = prepare_means_replicates(
                    piled_df, compartment_metabolites)
                dfs_dict = add_categorical_time(dfs_dict, time_levels_list)
                output_file_elements = [compartment, condition]

                stacked_plot_no_grid(
                    compartment_metabolites,
                    dfs_dict,
                    out_plot_dir,
                    output_file_elements,
                    cfg,
                    xlab_text=True,
                    x_to_plot="timenum",
                    x_ticks_text_tilt=x_ticks_text_tilt_fixed[
                                         'time_only'])
                plt.close()

                if (cfg.analysis.method.as_grid is not None) and \
                        cfg.analysis.method.as_grid:
                    complex_stacked_plot_as_grid(
                        compartment_metabolites, dfs_dict,
                        out_plot_dir,
                        output_file_elements,
                        cfg,
                        xlab_text=True,
                        x_to_plot="timenum",
                        x_ticks_text_tilt=x_ticks_text_tilt_fixed[
                            'time_only']
                    )

        elif (not cfg.analysis.method.split_plots_by_condition) and \
                (len(conditions) >= 2):
            if cfg.analysis.method.appearance_separated_time:
                conditions, combined_tc_levels = add_xemptyspace_tolabs(
                    conditions, time_levels_list)
            else:
                combined_tc_levels = time_plus_condi_labs(conditions,
                                                          time_levels_list)
            # end if
            dfs_dict = add_combined_conditime(dfs_dict, combined_tc_levels)
            output_file_elements = [compartment]

            stacked_plot_no_grid(compartment_metabolites,
                                 dfs_dict,
                                 out_plot_dir,
                                 output_file_elements,
                                 cfg,
                                 xlab_text=True,
                                 x_to_plot="time_and_condition",
                                 x_ticks_text_tilt=x_ticks_text_tilt_fixed[
                                     'time_and_condition'])
            plt.close()

            if (cfg.analysis.method.as_grid is not None) and \
                    cfg.analysis.method.as_grid:
                complex_stacked_plot_as_grid(
                    compartment_metabolites, dfs_dict, out_plot_dir,
                    output_file_elements,
                    cfg,
                    xlab_text=True,
                    x_to_plot="time_and_condition",
                    x_ticks_text_tilt=x_ticks_text_tilt_fixed[
                        'time_and_condition']
                )
            # end if
        # end if
        isotopologues_types_create_legend(cfg)  # legend alone
        plt.savefig(
            os.path.join(
                out_plot_dir, "legend_isotopologues_stackedbars.pdf"
            ),
            format="pdf")
        plt.close()
