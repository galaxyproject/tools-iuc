#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
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


def pile_up_abundance(df: pd.DataFrame,
                      metada_sel: pd.DataFrame) -> pd.DataFrame:
    dfcompartment = df.T.reset_index()
    dafull = pd.DataFrame(columns=["timepoint", "condition",
                                   "metabolite", "abundance"])

    for metabolite in dfcompartment.columns[1:]:
        subdf = dfcompartment[['index', metabolite]].rename(
            columns={'index': 'name_to_plot'})
        subdf = subdf.merge(metada_sel, on='name_to_plot')
        subdf['metabolite'] = metabolite
        subdf.rename(columns={metabolite: 'abundance'}, inplace=True)
        dafull = pd.concat([dafull, subdf[[
            'timepoint', 'condition', 'metabolite', 'abundance']]],
                           ignore_index=True)

    return dafull


def plot_one_metabolite(df: pd.DataFrame,
                        metabolite: str,
                        axisx_var: str,
                        hue_var: str,
                        axisx_labeltilt: int,
                        palette_choice: str,
                        curr_ax: matplotlib.axes
                        ) -> matplotlib.axes:
    """"
    returns a single object of type matplotlib.axes
    with all the individual metabolite plot
    """
    plt.rcParams.update({"font.size": 21})
    sns.barplot(
        ax=curr_ax,
        x=axisx_var,
        y="abundance",
        hue=str(hue_var),
        data=df,
        palette=palette_choice,
        alpha=1,
        edgecolor="black",
        errcolor="black",
        errwidth=1.7,
        capsize=0.12,
    )
    try:
        sns.stripplot(
            ax=curr_ax,
            x=axisx_var,
            y="abundance",
            hue=str(hue_var),
            data=df,
            palette=palette_choice,
            dodge=True,
            edgecolor="black",
            linewidth=1.5,
            alpha=1,
        )
    except Exception as e:
        # When Nan it throws error, avoid it
        print(e)
        pass

    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    curr_ax.set(title=metabolite + "\n")  # title of the subplot
    curr_ax.set(ylabel="Abundance")
    curr_ax.set(xlabel="")
    sns.despine()
    curr_ax.tick_params(axis="x", labelrotation=axisx_labeltilt)
    curr_ax.set_ylim(bottom=0)  # set minimal val display : y axis : 0

    return curr_ax


def plot_abundance_bars_no_grid(
    piled_sel_df: pd.DataFrame,
    selected_metabolites: List[str],
    CO: str,
    SMX: str,
    axisx_var: str,
    hue_var: str,
    output_directory: str,
    axisx_labeltilt: int,
    width_each_subfig: float,
    height_each_subfig: float,
    cfg: DictConfig,
) -> None:
    """ Saves independent single .pdf plot, one by each metabolite"""
    for k in range(len(selected_metabolites)):
        pile_df = piled_sel_df.loc[
                  piled_sel_df["metabolite"] == selected_metabolites[k], :]
        pile_df = pile_df.reset_index()
        output_path_k = os.path.join(
            output_directory,
            f"bars_{CO}_{selected_metabolites[k]}-{SMX}.pdf")

        fig_this_metabolite, axs_k = plt.subplots(
            nrows=1, ncols=1,
            figsize=(width_each_subfig, height_each_subfig))
        axs_k = plot_one_metabolite(pile_df, selected_metabolites[k],
                                    axisx_var, hue_var, axisx_labeltilt,
                                    cfg.analysis.method.palette,
                                    axs_k)
        if k == 0:  # collect the legend, which is same for any plot
            thehandles, thelabels = axs_k.get_legend_handles_labels()
        axs_k.legend_.remove()
        plt.tight_layout(pad=0.01, w_pad=-2, h_pad=0.1)
        plt.savefig(output_path_k, bbox_inches="tight", format="pdf")
        plt.close()

    plt.legend(handles=thehandles, labels=thelabels, loc="upper right")
    plt.axis("off")
    plt.savefig(os.path.join(output_directory, "legend.pdf"), format="pdf")
    plt.close()
    logger.info(f"Saved abundance plots in {output_directory}")


def plot_as_grid_of_bars(
        piled_sel_df: pd.DataFrame, selected_metabolites: List[str],
        CO: str, SMX: str,
        axisx_var: str, hue_var: str,
        output_directory: str,
        axisx_labeltilt: int,
        width_each_subfig: float,
        height_each_subfig: float,
        cfg: DictConfig
) -> None:
    """
    Saves as a grid (vertically arranged subplots)
    Note: corrector factors carefully tuned to match with size of independent
    plots in not grid version
    """
    corrector_factor = 0.7
    total_height_grid = height_each_subfig * len(selected_metabolites) + \
        (corrector_factor * 7 * len(selected_metabolites))

    fig, axs = plt.subplots(
        nrows=len(selected_metabolites), ncols=1,
        sharey=False, figsize=(width_each_subfig, total_height_grid))

    for i in range(len(selected_metabolites)):
        pile_df = piled_sel_df.loc[
                  piled_sel_df["metabolite"] == selected_metabolites[i], :]
        pile_df = pile_df.reset_index()

        axs[i] = plot_one_metabolite(pile_df, selected_metabolites[i],
                                     axisx_var, hue_var, axisx_labeltilt,
                                     cfg.analysis.method.palette,
                                     axs[i])

        if cfg.analysis.method.x_text_modify_as is not None:
            xticks_text_l: list = cfg.analysis.method.x_text_modify_as
            try:
                xticks_text_l = [str(i) for i in xticks_text_l]
                axs[i].set_xticklabels(xticks_text_l)
            except Exception as e:
                print(e, "The argument 'x_text_modify_as' is incorrectly set")

    thehandles, thelabels = axs[-1].get_legend_handles_labels()
    for i in range(len(selected_metabolites)):
        axs[i].legend_.remove()

    plt.subplots_adjust(left=0.2, top=0.76, bottom=0.2,
                        hspace=corrector_factor)

    output_path = os.path.join(output_directory, f"bars_{CO}_{SMX}.pdf")
    plt.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close()

    plt.legend(handles=thehandles, labels=thelabels, loc="upper right")
    plt.axis("off")
    plt.savefig(os.path.join(output_directory, "legend.pdf"), format="pdf")
    plt.close()
    logger.info(f"Saved abundance plot in {output_path}")


def run_plot_abundance_bars(dataset: Dataset, out_plot_dir,
                            cfg: DictConfig) -> None:
    metadata_df = dataset.metadata_df

    timepoints = cfg.analysis.timepoints  # locate where it is used
    metabolites = (
        cfg.analysis.metabolites
    )  # will define which metabolites are plotted in the abundance plot
    conditions = cfg.analysis.dataset.conditions  # <= locate where it is used

    axisx_labeltilt = cfg.analysis.method.axisx_labeltilt
    axisx_var = cfg.analysis.method.axisx
    hue_var = cfg.analysis.method.barcolor

    width_each_subfig = cfg.analysis.width_each_subfig
    height_each_subfig = cfg.analysis.method.height_each_subfig

    compartments = set(metadata_df["short_comp"])
    for compartment in compartments:
        metadata_compartment_df: pd.DataFrame = \
            metadata_df.loc[metadata_df["short_comp"] == compartment, :]
        compartment_df = dataset.compartmentalized_dfs[
            "abundances"][compartment]
        # metadata and abundances time of interest
        metadata_slice = metadata_compartment_df.loc[
                         metadata_compartment_df[
                             "timepoint"].isin(timepoints), :]
        values_slice = compartment_df[list(metadata_slice["name_to_plot"])]

        # total piled-up data:
        piled_sel = pile_up_abundance(values_slice, metadata_slice)
        piled_sel["condition"] = pd.Categorical(
            piled_sel["condition"], conditions)
        piled_sel["timepoint"] = pd.Categorical(
            piled_sel["timepoint"], timepoints)

        plot_abundance_bars_no_grid(
            piled_sel,
            metabolites[compartment],
            compartment,
            "total_abundance",
            axisx_var,
            hue_var,
            out_plot_dir,
            axisx_labeltilt,
            width_each_subfig,
            height_each_subfig,
            cfg
        )

        if cfg.analysis.method.as_grid is not None:
            if cfg.analysis.method.as_grid:
                plot_as_grid_of_bars(
                    piled_sel, metabolites[compartment], compartment,
                    "total_abundance",
                    axisx_var, hue_var,
                    out_plot_dir, axisx_labeltilt, width_each_subfig,
                    height_each_subfig, cfg
                )
