#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
from typing import List, Union

from hydra.core.config_store import ConfigStore

import matplotlib.figure as figure
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import numpy as np

from omegaconf import DictConfig

import pandas as pd

import seaborn as sns


logger = logging.getLogger(__name__)

cs = ConfigStore.instance()


def variance_expl_plot(var_explained_df: pd.DataFrame) -> figure.Figure:
    """
    returns the bar-plot with the percentage of explained variances by PC
    """
    fig = plt.figure()
    sns.barplot(x='PC', y='Explained Variance %',
                data=var_explained_df, color="cadetblue")
    plt.title("Percent variability explained by the principal components")
    plt.ylabel("Explained variances (%)")
    return fig


def eigsorted(cov: np.array) -> tuple:
    """
    used for calculating ellipses
    many thanks to :
    https://rayblick.gitbooks.io/my-python-scrapbook/content/
    analysis/plotting/scatterplot_ellipse.html
    """
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def pca_scatter_plot(pc_df: pd.DataFrame,
                     var_explained_df: pd.DataFrame,
                     col1: str, col2: str,
                     labels_column: str,
                     ellipses_column: Union[str, None]) -> figure.Figure:
    """
    returns the scatter plot (a seaborn matplotlib figure)
    """
    fig, ax = plt.subplots()
    sns.scatterplot(x="PC1", y="PC2",
                    ax=ax,
                    data=pc_df,
                    hue=col1,
                    style=col2,
                    legend=True,
                    s=80, zorder=3)
    ax.axhline(0, ls="--", color="gray", zorder=1)
    ax.axvline(0, ls="--", color="gray", zorder=1)
    if labels_column != "":
        for i, row in pc_df.iterrows():
            ax.text(pc_df.at[i, 'PC1'] + 0.2, pc_df.at[i, 'PC2'],
                    pc_df.at[i, labels_column], size='x-small')
    # end if
    row_xlab = var_explained_df.iloc[0, :]
    row_ylab = var_explained_df.iloc[1, :]

    plt.xlabel(
        f"{row_xlab['PC']} {round(row_xlab['Explained Variance %'], 2)} %")
    plt.ylabel(
        f"{row_ylab['PC']} {round(row_ylab['Explained Variance %'], 2)} %")
    plt.title("")

    if ellipses_column is not None:   # add ellipses if specified in *args
        ellipses_groups_list = pc_df[ellipses_column].unique()
        for group in ellipses_groups_list:
            xdata = pc_df.loc[pc_df[ellipses_column] == group, 'PC1']
            ydata = pc_df.loc[pc_df[ellipses_column] == group, 'PC2']
            # get values to build the ellipse
            cov = np.cov(xdata, ydata)
            vals, vecs = eigsorted(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            w, h = 2 * 2 * np.sqrt(vals)
            # create the ellipse
            ell = Ellipse(xy=(np.mean(xdata), np.mean(ydata)),
                          width=w, height=h,
                          angle=theta,
                          edgecolor="lightgray",
                          linestyle='-', facecolor='none')
            ax.add_artist(ell)

    return fig


def pca_scatter_2_pdf(figure_pc: figure.Figure,
                      name_elements: List[str],
                      out_plot_dir: str) -> None:
    name_plot = f"{'--'.join(name_elements)}_pc.pdf"
    figure_pc.savefig(os.path.join(out_plot_dir, name_plot))


def demo_pca_iris(out_plot_dir: str) -> None:
    from sklearn.decomposition import PCA
    iris = sns.load_dataset("iris")
    sns.relplot(data=iris, x="sepal_width", y="petal_width",
                hue="species")
    iris = iris.assign(name_to_plot=[str(i) for i in iris.index])
    iris_metadata = iris[['name_to_plot', "species"]]
    iris = iris.drop(columns=['name_to_plot', "species"])
    df = iris.T  # variables rows, samples columns
    df = df.div(df.std(axis=1, ddof=0), axis=0)
    df.columns = [str(i) for i in iris.index]
    X = np.transpose(np.array(df))
    pca = PCA(n_components=4)
    pc = pca.fit_transform(X)
    pc_df = pd.DataFrame(data=pc,
                         columns=['PC' + str(i) for i in range(1, 4 + 1)])
    pc_df = pc_df.assign(name_to_plot=df.columns)
    pc_df = pd.merge(pc_df, iris_metadata, on='name_to_plot')
    var_explained_df = pd.DataFrame({
        'Explained Variance %': pca.explained_variance_ratio_ * 100,
        'PC': ['PC' + str(i) for i in range(1, 4 + 1)]})

    name_elements = ['IrisDemo']
    scatter_fig = pca_scatter_plot(pc_df,
                                   var_explained_df, "species",
                                   "species", labels_column="",
                                   ellipses_column="species")

    pca_scatter_2_pdf(scatter_fig, name_elements, out_plot_dir)
    logger.info("Saving demo figure: pca on iris dataset")
    plt.close()


def run_pca_plot(pca_results_dict: dict,  cfg: DictConfig,
                 out_plot_dir: str) -> None:
    for tup in pca_results_dict.keys():
        pc_df = pca_results_dict[tup]['pc']
        var_explained_df = pca_results_dict[tup]['var']
        figure_var: figure.Figure = variance_expl_plot(var_explained_df)
        name_plot_var = f"{'--'.join(tup)}_var.pdf"
        figure_var.savefig(os.path.join(out_plot_dir, name_plot_var))
        plt.close()
        options_labels = {'label-y': "name_to_plot",
                          'label-n': ""}  # when empty string, no dot labels
        # scatter: save both versions, labeled dots and unlabeled dots:
        for choice in options_labels.keys():
            labels_column = options_labels[choice]
            name_elements = list(tup) + [choice]
            scatter_fig: figure.Figure = pca_scatter_plot(
                pc_df,  var_explained_df, "condition",
                "condition", labels_column,
                ellipses_column=cfg.analysis.method.draw_ellipses)
            pca_scatter_2_pdf(scatter_fig, name_elements, out_plot_dir)
            plt.close()

    # end for
    if cfg.analysis.method.run_iris_demo:
        demo_pca_iris(out_plot_dir)

# end

# calculate Ellipses , many thanks to :
# https://rayblick.gitbooks.io/my-python-scrapbook/content/analysis/plotting/scatterplot_ellipse.html

# ex4:
# https://www.programcreek.com/python/example/61396/matplotlib.patches.Ellipse
