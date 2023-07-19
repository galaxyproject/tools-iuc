#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski

contexts: mean the comparisons in this version, but evolution should
   handle concentrations/expressions/..., not only the differential results
"""
import logging
import os
import warnings
from typing import Dict, List

from constants import (
    assert_literal,
    availtest_methods_type,
    data_files_keys_type,
    molecular_types_for_metabologram
)

from data import DataIntegration

import helpers

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize

import numpy as np

from omegaconf import DictConfig

import pandas as pd

from processing import differential_analysis

import seaborn as sns


logger = logging.getLogger(__name__)


def get_differential_results_dict(file_name: str,
                                  data_integration: DataIntegration,
                                  cfg: DictConfig,
                                  test: str,
                                  compartment: str
                                  ) -> Dict[int, pd.DataFrame]:
    impute_value = cfg.analysis.method.impute_values[file_name]
    df = data_integration.compartmentalized_dfs[file_name][compartment]

    df = df[(df.T != 0).any()]
    val_instead_zero = helpers.arg_repl_zero2value(impute_value, df)
    df = df.replace(to_replace=0, value=val_instead_zero)

    metabolomics_differential_results_dict = {}  # Dict[int, pd.DataFrame]

    for i, comparison in enumerate(cfg.analysis.comparisons):
        result = differential_analysis.pairwise_comparison(
            df, data_integration, cfg, comparison, test)
        result["compartment"] = compartment
        result = differential_analysis.reorder_columns_diff_end(result)
        result = result.sort_values(["padj", "distance/span"],
                                    ascending=[True, False])
        result.reset_index(names=[cfg.analysis.columns_metabolites['ID']],
                           inplace=True)
        metabolomics_differential_results_dict[i] = result

    return metabolomics_differential_results_dict


def pathways_df_2_lists_dict(pathways_df: pd.DataFrame
                             ) -> Dict[str, List[str]]:
    """
    transform dictionary of dataframes,
    into nested a dictionary of lists like this:
    {
    metabolites:{'PATH1': ['Citric_acid', 'Fumaric_acid'], 'PATH2':...},
    transcripts: {'PATH1': ['FOX', 'ILE5', ...], 'PATH2': ....}
    }
    """
    pathways_df = pathways_df.fillna("Empty")
    preD = pathways_df.to_dict(orient='list')
    pathD = dict()
    for k in preD.keys():
        deduplicated = set(list(preD[k])) - set(["Empty", "", np.nan])
        pathD[k] = list(deduplicated)
    return pathD


def pile_dfs_by_contexts(contexts_dict: Dict[int, Dict[str, pd.DataFrame]],
                         columns_dict: Dict[str, Dict[str, str]],
                         ) -> Dict[str, pd.DataFrame]:
    """
    returns dictionary of dataframes:
    - key is the type of molecule (e.g. metabolites),
    - value is the dataframe:
     having piled up the contexts (e.g. comparisons) in column context:
            name      VALUES   context   typemol
            Citrate    -1.11      0         metabolites
            Fumarate   4.08       1         metabolites
            ...
    """
    data_cleaned_dict = {}
    for molecule_type in molecular_types_for_metabologram:
        df_out = pd.DataFrame({'name': [], 'VALUES': []})
        for context in contexts_dict.keys():
            df = contexts_dict[context][molecule_type]
            df = df[[columns_dict[molecule_type]['ID'],
                     columns_dict[molecule_type]['values']]]
            df.columns = ['name', 'VALUES']
            df = df.assign(VALUES=df['VALUES'].round(4))
            df = df.assign(context=context, typemol=molecule_type)
            df['context'] = df['context'].astype(int)
            df_out = pd.concat([df_out, df], axis=0)

        data_cleaned_dict[molecule_type] = df_out
    return data_cleaned_dict


def filter_by_pathway_dict(data_cleaned_dict, pathways_lists_dict):
    """
    molecules that are both in pathways and data are kept
    """
    for type_molecule in data_cleaned_dict.keys():
        df = data_cleaned_dict[type_molecule].copy()
        pathways_molecules = pathways_lists_dict[type_molecule].values()
        aggregated_pathways_by_type_molecule = set()
        for list_of_molecules in pathways_molecules:
            for i in list_of_molecules:
                aggregated_pathways_by_type_molecule.update([i])

        df = df.loc[df.name.isin(
            list(aggregated_pathways_by_type_molecule)), :]
        data_cleaned_dict[type_molecule] = df

    return data_cleaned_dict


def set_absolute_max_by_type_molecule(data_cleaned_dict, cfg
                                      ) -> Dict[str, float]:
    max_absolute_value_dict = {}
    for k in cfg.analysis.method.abs_values_scale_color_bar.keys():
        # keys : metabolites, transcripts
        if cfg.analysis.method.abs_values_scale_color_bar[k] is None:
            max_absolute_value_dict[k] = max(
                abs(data_cleaned_dict[k]['VALUES']))
        else:
            try:
                v = float(cfg.analysis.method.abs_values_scale_color_bar[k])
            except ValueError:
                v = max(abs(data_cleaned_dict[k]['VALUES']))
            max_absolute_value_dict[k] = v
    return max_absolute_value_dict


def get_custom_color_palette_hash(lowcolor: str, midcolor: str,
                                  highcolor: str):
    """
    courtesy from :
    https://richardhildebrand.wordpress.com/2019/09/18/create-a-custom-color-palette-with-matplotlib-and-seaborn/
    """
    colorlist = [lowcolor, midcolor, highcolor]
    return LinearSegmentedColormap.from_list("", colorlist, N=256)


def rgbas2hex(rgbas_list: List[tuple]):
    colorsout = []
    for tup in rgbas_list:
        tmp = matplotlib.colors.to_hex(tup)
        colorsout.append(tmp)
    return colorsout


def values2rgbas(myvalues: np.array, my_cmap: LinearSegmentedColormap,
                 vmin: float, vmax: float, center: float) -> List[tuple]:
    if center == 0:
        # Normalize data before giving colors,
        # because map interval is [0,1] by matplotlib
        # https://stackoverflow.com/questions/25408393/getting-individual-colors-from-a-color-map-in-matplotlib
        norm = Normalize(vmin=vmin, vmax=vmax)
        rgba_tuples = my_cmap(norm(myvalues))
        return rgba_tuples
    else:
        print("only center == 0 is handled here")


def set_colormap(data_cleaned_dict, my_cmap, max_absolute_value_dict):
    for molecule_type in data_cleaned_dict.keys():
        v = max_absolute_value_dict[molecule_type]
        data_cleaned_dict[molecule_type]["mycolors"] = rgbas2hex(
            values2rgbas(
                data_cleaned_dict[molecule_type]['VALUES'].to_numpy(),
                my_cmap, -v, v, center=0)
        )
    return data_cleaned_dict


def cleaned_dict_2_df(data_cleaned_dict, cfg
                      ) -> (pd.DataFrame, Dict[str, float], float):
    """
    returns dataframe, piled transcripts and metabolites data:
    name      VALUES   context   typemol      mycolors
    Citrate    -1.11      0         metabolites   #ffa500
    ABCA1      4.08       1         metabolites   #80c5c5
    ...
    """
    max_absolute_value_dict = set_absolute_max_by_type_molecule(
        data_cleaned_dict, cfg)
    my_cmap = get_custom_color_palette_hash(
        lowcolor=cfg.analysis.method.colors_divergent_palette[0],
        midcolor=cfg.analysis.method.colors_divergent_palette[1],
        highcolor=cfg.analysis.method.colors_divergent_palette[2])
    # color scale differs metabolites or transcripts
    data_cleaned_dict = set_colormap(data_cleaned_dict, my_cmap,
                                     max_absolute_value_dict)

    # gathered_df
    gathered = pd.concat([data_cleaned_dict["metabolites"],
                          data_cleaned_dict["transcripts"]], axis=0)

    gathered["molecule_label"] = gathered["name"]
    values_str_list = [str(i) for i in gathered["VALUES"]]

    if cfg.analysis.method.display_label_and_value:
        gathered["molecule_label"] = gathered["name"].str.cat(values_str_list,
                                                              sep=": ")

    return gathered, max_absolute_value_dict, my_cmap


def get_ordered__metabolome_contexts(cfg: DictConfig) -> Dict[int, str]:
    ordered__metabolome_contexts: Dict[int, str] = {}
    for i, context in enumerate(cfg.analysis.comparisons):
        tmp = "-".join("-".join([ele for ele in sublist])
                       for sublist in context)
        ordered__metabolome_contexts[i] = tmp

    return ordered__metabolome_contexts


def metabologram_organize_data(
        dam_dfs_dict: Dict[int, pd.DataFrame],
        deg_dfs_dict: Dict[int, pd.DataFrame],
        titles_dict: Dict[int, str],
        pathways_lists_dict: Dict[str, Dict[str, List[str]]],
        cfg: DictConfig
) -> (
        pd.DataFrame, dict, LinearSegmentedColormap
):
    """refactors and cleans data dictionaries"""

    columns_dict = {'metabolites': {
        'ID': cfg.analysis.columns_metabolites['ID'],
        'values': cfg.analysis.columns_metabolites['values']
    },
        'transcripts': {
            'ID': cfg.analysis.columns_transcripts['ID'],
            'values': cfg.analysis.columns_transcripts['values']
        }
    }

    contexts_dict = dict()
    for k in titles_dict.keys():  # k is integer
        contexts_dict[k] = {'transcripts': deg_dfs_dict[k],
                            'metabolites': dam_dfs_dict[k],
                            'title': titles_dict[k]}

    data_cleaned_dict = pile_dfs_by_contexts(contexts_dict, columns_dict)
    data_cleaned_dict = filter_by_pathway_dict(
        data_cleaned_dict, pathways_lists_dict)
    gathered, max_absolute_value_dict, my_cmap = cleaned_dict_2_df(
        data_cleaned_dict, cfg)

    return gathered, max_absolute_value_dict, my_cmap


# # viz functions


def inner_pie_colors(inner_dict, my_cmap,
                     max_abs_transcripts: float, max_abs_metabolites: float):
    metaboval = inner_dict['metabo_mean_val']
    geneval = inner_dict['gene_mean_val']
    metabocolors_l = rgbas2hex(
        values2rgbas(
            [metaboval], my_cmap, -max_abs_metabolites, max_abs_metabolites,
            center=0)
    )
    genecolors_l = rgbas2hex(
        values2rgbas(
            [geneval], my_cmap, -max_abs_transcripts, max_abs_transcripts,
            center=0)
    )
    return {'metab_color': metabocolors_l[0], 'gene_color': genecolors_l[0]}


def introduce_nan_elems_if_not_in(curr_pathway_df: pd.DataFrame,
                                  path_elems_here: list,
                                  context_here: str,
                                  genes_list: list,
                                  metabo_list: list) -> pd.DataFrame:
    """"
    for the current pathway, fill with nan in df where no matches
    """
    not_in_data = set(path_elems_here) - set(curr_pathway_df['name'].tolist())
    nan_df = pd.DataFrame(columns=['name', 'VALUES', 'context', 'typemol',
                                   'mycolors', 'molecule_label'])
    nan_df = nan_df.assign(name=list(not_in_data))
    nan_df = nan_df.assign(VALUES=np.nan)
    nan_df = nan_df.assign(context=context_here)
    nan_df = nan_df.assign(typemol='')
    nan_df = nan_df.assign(mycolors="gray")
    nan_df = nan_df.assign(molecule_label=list(not_in_data))

    for i, r in nan_df.iterrows():
        if r['name'] in genes_list:
            nan_df.loc[i, 'typemol'] = 'transcripts'
        elif r['name'] in metabo_list:
            nan_df.loc[i, 'typemol'] = 'metabolites'

    curr_pathway_df = pd.concat([curr_pathway_df, nan_df], axis=0)

    return curr_pathway_df


def donut_outer(curr_pathway_context_df,
                cfg: DictConfig, fig: matplotlib.figure.Figure
                ) -> matplotlib.figure.Figure:
    """ external portion of the donut plot """
    curr_pathway_context_df['circportion'] = ''
    genecircportion = 50 / curr_pathway_context_df.loc[
                           curr_pathway_context_df.typemol == "transcripts",
                           :].shape[0]
    metabocircportion = 50 / curr_pathway_context_df.loc[
                             curr_pathway_context_df.typemol == "metabolites",
                             :].shape[0]
    curr_pathway_context_df.loc[
        curr_pathway_context_df.typemol == "transcripts",
        "circportion"] = genecircportion
    curr_pathway_context_df.loc[
        curr_pathway_context_df.typemol == "metabolites",
        "circportion"] = metabocircportion

    sizes_list = curr_pathway_context_df["circportion"]
    annots = curr_pathway_context_df["molecule_label"]
    mappedcolors_list = curr_pathway_context_df["mycolors"]

    plt.pie(sizes_list,
            colors=mappedcolors_list,
            wedgeprops={'width': 1,
                        'edgecolor': cfg.analysis.method.edge_color[0],
                        'linewidth': cfg.analysis.method.line_width[0]},
            radius=1,
            startangle=90,
            labels=annots,
            # this one yiels the  labels annotated in the plot
            textprops={'fontsize': cfg.analysis.method.font_size})
    # white circles for artist patches
    ax = fig.add_subplot()

    ax.add_patch(plt.Circle((0, 0), radius=0.49,
                            edgecolor=cfg.analysis.method.edge_color[0],
                            linewidth=1.6))
    ax.add_patch(plt.Circle((0, 0), radius=0.482, color="white"))
    return fig


def donut_inner(gatheredsub, cfg: DictConfig,
                my_cmap, max_absolute_value_dict: Dict[str, float],
                fig: matplotlib.figure.Figure) -> matplotlib.figure.Figure:
    """central part of the donut plot"""
    inner_dict = {'metabo_mean_val': gatheredsub.loc[
        gatheredsub.typemol == 'metabolites', 'VALUES'].mean(),
                  'gene_mean_val': gatheredsub.loc[
                      gatheredsub.typemol == 'transcripts', 'VALUES'].mean()}
    inner_colorsD = inner_pie_colors(inner_dict, my_cmap,
                                     max_absolute_value_dict["transcripts"],
                                     max_absolute_value_dict["metabolites"])
    # internal pie
    plt.pie([50, 50],
            colors=[inner_colorsD['metab_color'],
                    inner_colorsD['gene_color']],
            wedgeprops={'width': 0.41,
                        'edgecolor': cfg.analysis.method.edge_color[1],
                        'linewidth': cfg.analysis.method.line_width[1]},
            radius=0.41,
            startangle=90,
            labels=np.array([inner_dict['metabo_mean_val'].round(2),
                             inner_dict['gene_mean_val']]).round(1),
            labeldistance=0.2)
    return fig


def combine_elements_by_pathway_2dict(
        pathways_lists_dict: Dict[str, Dict[str, List[str]]]
) -> Dict[str, List[str]]:
    """
    returns dictionary :
    {'PATH1': [metabolite1, metabolite2, .., GENE1, GENE2, ...],
    'PATH2': ['metabolite3', 'Metabolite4', .., GENE3, GENE4, ...]}
    """
    combined_elements_by_pathway = {}  # finalD
    for type_of_data in pathways_lists_dict.keys():
        for pathway in pathways_lists_dict[type_of_data].keys():
            if pathway not in combined_elements_by_pathway.keys():
                combined_elements_by_pathway[pathway] = \
                    pathways_lists_dict[type_of_data][pathway].copy()
            else:
                combined_elements_by_pathway[pathway] += (
                    [i for i in pathways_lists_dict[type_of_data][pathway]]
                )
    return combined_elements_by_pathway


def save_donuts_plots(gathered: pd.DataFrame,
                      max_absolute_value_dict,
                      my_cmap, cfg: DictConfig,
                      pathways_lists_dict: dict,
                      titles_dict: dict, compartment: str, file_name: str,
                      out_plot_dir: str
                      ) -> None:
    combined_elements_by_pathway = combine_elements_by_pathway_2dict(
        pathways_lists_dict
    )

    # prepare auxiliary_indexes_figs_dict to avoid doing nested loops
    auxiliary_indexes_figs_dict = dict()
    indexer = 0
    for i in combined_elements_by_pathway.keys():
        for j in titles_dict.keys():
            auxiliary_indexes_figs_dict[indexer] = {
                'path': i,
                'context': j,
                'title': [str(i), titles_dict[j], file_name]}
            indexer += 1

    for indexer, auxiliary_inner_dict \
            in enumerate(auxiliary_indexes_figs_dict.values()):
        pathway_k: str = auxiliary_inner_dict["path"]
        path_elems_here = combined_elements_by_pathway[pathway_k]
        gathered_sub = gathered.loc[gathered['name'].isin(path_elems_here), :]
        context_here = auxiliary_indexes_figs_dict[indexer]['context']
        gathered_sub = introduce_nan_elems_if_not_in(
            gathered_sub,
            path_elems_here,
            context_here,
            pathways_lists_dict["transcripts"][pathway_k],
            pathways_lists_dict["metabolites"][pathway_k])

        gathered_sub['typemol'] = pd.Categorical(gathered_sub['typemol'],
                                                 categories=['metabolites',
                                                             'transcripts'])
        gathered_sub = gathered_sub.sort_values(by=['typemol', 'name'],
                                                ascending=[True, False])
        gathered_sub = gathered_sub.loc[
                       gathered_sub['context'] == context_here, :]

        fig = plt.figure(figsize=(cfg.analysis.method.fig_width,
                                  cfg.analysis.method.fig_height))
        fig.tight_layout()
        fig.suptitle(f"{auxiliary_inner_dict['title'][1]}\n"
                     f"{auxiliary_inner_dict['title'][0]}")

        fig = donut_outer(gathered_sub, cfg, fig)
        fig = donut_inner(gathered_sub, cfg, my_cmap,
                          max_absolute_value_dict, fig)
        # end donut
        title_out = "-".join(auxiliary_inner_dict['title'])
        base_file_name = f"{title_out}-{compartment}"
        out_path = os.path.join(out_plot_dir, base_file_name)
        plt.savefig(f"{out_path}.pdf")

    # end for

    # Legends: color key bars
    fig, axes = plt.subplots(ncols=2, nrows=1,
                             figsize=(cfg.analysis.method.fig_width,
                                      cfg.analysis.method.fig_height))
    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore")  # Remove warning internal heatmap (df.shape)
        vs = (max_absolute_value_dict["metabolites"],
              max_absolute_value_dict["transcripts"])
        labels = ('Metabolite', 'Transcript')
        for ax, v, label in zip(axes, vs, labels):
            sns.heatmap([[]], ax=ax, cmap=my_cmap, center=0, cbar=True,
                        annot=False, yticklabels=False,
                        square=True,
                        vmin=-v, vmax=v,
                        cbar_kws={'shrink': 0.9, 'aspect': 10,
                                  'label': label,
                                  'drawedges': False})

        base_file_name = f"legend-{file_name}-{compartment}"
        out_path = os.path.join(out_plot_dir, base_file_name)
        plt.savefig(f"{out_path}.pdf")


def run_metabologram(file_name: str,
                     data_integration: DataIntegration,
                     cfg: DictConfig,
                     test: str, out_plot_dir: str):
    assert_literal(test, availtest_methods_type, "Available test")
    assert_literal(file_name, data_files_keys_type, "file name")
    compartment: str = cfg.analysis.compartment  # one compartment admitted
    deg_dfs_dict: Dict[int, pd.DataFrame] = data_integration.deg_dfs
    pathways_dfs_dict: Dict[str, pd.DataFrame] = data_integration.pathways_dfs
    ordered__metabolome_contexts = get_ordered__metabolome_contexts(cfg)
    names_transcripts_dict = data_integration.get_names_transcripts_files()

    titles_dict: Dict[int, str] = {}
    for i in names_transcripts_dict.keys():
        # keys are simply integers (user order) for the following dicts:
        titles_dict[i] = "--".join([ordered__metabolome_contexts[i],
                                    names_transcripts_dict[i]])

    logger.info(f"Metabologram running for the user specified integration of"
                f" metabolomics-transcriptomics : {titles_dict}")

    dam_dfs_dict: Dict[int, pd.DataFrame] = \
        get_differential_results_dict(file_name, data_integration,
                                      cfg, test, compartment)

    ############
    # organize pathways
    pathways_lists_dict = {}
    for type_of_data in pathways_dfs_dict.keys():
        # keys are : transcripts, metabolites
        pathways_lists_dict[type_of_data] = pathways_df_2_lists_dict(
            pathways_dfs_dict[type_of_data])
    ############
    # organize data
    gathered, max_absolute_value_dict, my_cmap = metabologram_organize_data(
        dam_dfs_dict, deg_dfs_dict, titles_dict,
        pathways_lists_dict, cfg)

    #########
    # plot
    save_donuts_plots(gathered, max_absolute_value_dict, my_cmap, cfg,
                      pathways_lists_dict,
                      titles_dict, compartment, file_name, out_plot_dir)

    logger.info(f"Saved the results in {out_plot_dir}")
