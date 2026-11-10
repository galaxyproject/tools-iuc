#!/usr/bin/env python
"""Run downstream TISSUE analyses on calibrated spatial AnnData."""

import argparse
from importlib.metadata import version

import anndata
import numpy as np
import pandas as pd
import tissue.downstream


def interval_widths(data, predicted, normalize):
    widths = data.obsm[predicted + "_hi"].values - data.obsm[predicted + "_lo"].values
    if not np.isfinite(widths).all() or (widths < 0).any():
        raise ValueError("Prediction intervals must have finite, nonnegative widths for filtering or PCA.")
    if normalize:
        widths = widths / (1 + np.abs(data.obsm[predicted].values))
    return widths


def obs_values(data, key):
    if key not in data.obs:
        raise ValueError(f"Observation metadata column {key!r} is absent from the input.")
    if data.obs[key].isna().any():
        raise ValueError(f"Observation metadata column {key!r} contains missing values.")
    return data.obs[key].astype(str).values


def run(args):
    data = anndata.read_h5ad(args.input)
    predicted = f"{args.method}_predicted_expression"
    if predicted not in data.obsm:
        raise ValueError(f"Input must contain calibrated predictions in obsm[{predicted!r}].")
    if not np.isfinite(data.obsm[predicted].values).all():
        raise ValueError("Predicted expression must be finite.")
    np.random.seed(args.seed)

    if args.operation == "testing":
        for suffix in ("_score", "_error", "_uncertainty", "_groups"):
            if predicted + suffix not in data.obsm:
                raise ValueError(f"Input is missing calibration metadata: {predicted + suffix}.")
        if "conf_genes_used" not in data.uns:
            raise ValueError("Input is missing uns['conf_genes_used'].")
        if args.test == "spatialde":
            if "spatial" not in data.obsm:
                raise ValueError("The SpatialDE test requires spatial coordinates in obsm['spatial'].")
            condition = None
        else:
            groups = obs_values(data, args.condition)
            data.obs[args.condition] = groups
            if args.group2 and not args.group1:
                raise ValueError("A second comparison group requires a first comparison group.")
            first_groups = [args.group1] if args.group1 else np.unique(groups)
            for first in first_groups:
                second_mask = groups == args.group2 if args.group2 else groups != first
                if (groups == first).sum() < 2 or second_mask.sum() < 2:
                    raise ValueError("Each comparison group must contain at least two cells.")
                if args.group2 == first:
                    raise ValueError("Comparison groups must be distinct.")
            condition = args.condition
        keys = tissue.downstream.multiple_imputation_testing(
            data, predicted, calib_genes=list(data.uns["conf_genes_used"]),
            condition=condition, test=args.test, n_imputations=args.n_imputations,
            group1=args.group1, group2=args.group2, symmetric=args.symmetric, return_keys=True,
        )
        tables = []
        for key in keys:
            result = data.uns[key]
            tables.append(pd.DataFrame({"result": key, "gene": result.columns, "value": result.iloc[0].values}))
        pd.concat(tables, ignore_index=True).to_csv(args.table, sep="\t", index=False, float_format="%.10g")
    elif args.operation == "weighted_pca":
        interval_widths(data, predicted, False)
        if args.n_components > min(data.obsm[predicted].shape):
            raise ValueError("PCA components cannot exceed the predicted gene count or cell count.")
        binarize = {"none": False, "otsu": True, "quantile": args.binarize_quantile}[args.binarize]
        weights = tissue.downstream.weighted_PCA(
            data, args.method, pca_method=args.pca_method, weighting=args.weighting,
            quantile_cutoff=args.quantile_cutoff, n_components=args.n_components,
            replace_inf=None if args.replace_inf == "none" else args.replace_inf,
            binarize=binarize, binarize_ratio=args.binarize_ratio,
            log_transform=args.log_transform, scale=args.scale, return_weights=True,
        )
        data.obsm["tissue_pca_weights"] = weights
        write_pcs(data, predicted, args)
    else:
        widths = interval_widths(data, predicted, args.normalize)
        strata = obs_values(data, args.stratification) if args.stratification else None
        proportion = "otsu" if args.filter_method == "otsu" else args.proportion
        # Constant gene-wise widths give undefined z-scores in the upstream filter.
        if strata is None:
            subsets = [widths]
        else:
            subsets = [widths[strata == label] for label in np.unique(strata)]
        if any((np.std(subset, axis=0) == 0).any() for subset in subsets):
            raise ValueError("Filtering requires variable interval widths for every gene within each stratum.")
        keep = tissue.downstream.detect_uncertain_cells(widths, proportion=proportion, stratification=strata)
        if args.operation == "filtered_pca":
            if args.n_components > min(len(keep), data.obsm[predicted].shape[1]):
                raise ValueError("PCA components cannot exceed the retained cell count or predicted gene count.")
            keep = tissue.downstream.filtered_PCA(
                data, args.method, proportion=proportion, stratification=strata,
                n_components=args.n_components, scale=args.scale, normalize=args.normalize,
                return_keep_idxs=True,
            )
        data.obs["tissue_keep"] = np.isin(np.arange(data.n_obs), keep)
        if args.operation == "filtered_pca":
            write_pcs(data, predicted, args)
        else:
            data.obs[["tissue_keep"]].to_csv(args.table, sep="\t", index_label="cell_id")

    data.uns["tissue_downstream_galaxy"] = {
        "tissue_version": version("tissue-sc"), "operation": args.operation,
        "method": args.method, "seed": args.seed,
    }
    data.write_h5ad(args.output, compression="gzip")


def write_pcs(data, predicted, args):
    pcs = data.obsm[f"{predicted}_PC{args.n_components}_"]
    if not np.isfinite(pcs).all():
        raise ValueError("PCA produced nonfinite components; check interval widths and weighting settings.")
    frame = pd.DataFrame(pcs, index=data.obs_names, columns=[f"PC{i + 1}" for i in range(args.n_components)])
    if args.operation == "filtered_pca":
        frame.insert(0, "tissue_keep", data.obs["tissue_keep"])
    frame.to_csv(args.table, sep="\t", index_label="cell_id", float_format="%.10g")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=version("tissue-sc"))
    parser.add_argument("--input", required=True)
    parser.add_argument("--method", choices=["spage", "knn", "tangram"], default="spage")
    parser.add_argument("--operation", choices=["testing", "weighted_pca", "filtered_pca", "detect"], required=True)
    parser.add_argument("--seed", type=int, default=444)
    parser.add_argument("--condition")
    parser.add_argument("--group1")
    parser.add_argument("--group2")
    parser.add_argument("--test", choices=["ttest", "wilcoxon_less", "wilcoxon_greater", "spatialde"], default="ttest")
    parser.add_argument("--n-imputations", type=int, default=100)
    parser.add_argument("--symmetric", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--n-components", type=int, default=15)
    parser.add_argument("--scale", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--pca-method", choices=["wpca", "empca", "pca"], default="wpca")
    parser.add_argument("--weighting", choices=["uniform", "inverse_pi_width", "inverse_norm_pi_width"], default="inverse_norm_pi_width")
    parser.add_argument("--quantile-cutoff", type=float)
    parser.add_argument("--replace-inf", choices=["none", "max", "min", "mean", "median"], default="none")
    parser.add_argument("--binarize", choices=["none", "otsu", "quantile"], default="quantile")
    parser.add_argument("--binarize-quantile", type=float, default=0.2)
    parser.add_argument("--binarize-ratio", type=float, default=10)
    parser.add_argument("--log-transform", action="store_true")
    parser.add_argument("--filter-method", choices=["proportion", "otsu"], default="proportion")
    parser.add_argument("--proportion", type=float, default=0.05)
    parser.add_argument("--stratification")
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args()
    if args.operation == "testing" and args.test != "spatialde" and not args.condition:
        parser.error("--condition is required for hypothesis testing")
    if args.operation == "testing" and args.test == "spatialde" and (args.group1 or args.group2):
        parser.error("comparison groups are not supported with the SpatialDE test")
    try:
        run(args)
    except (ValueError, KeyError, np.linalg.LinAlgError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
