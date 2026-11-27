#!/usr/bin/env python
"""Run the TISSUE Python API and export its results for Galaxy."""

import argparse
from importlib.metadata import version

import anndata
import numpy as np
import tissue.main
from scipy import sparse

__version__ = "1.0.1"


def read_data(path, label):
    # Normalization and upstream prediction/calibration require in-memory X.
    data = anndata.read_h5ad(path)
    # Match genes before selecting shared genes and withholding targets.
    data.var_names = data.var_names.str.lower()
    if not data.var_names.is_unique:
        raise ValueError(f"{label} gene names must be unique after conversion to lowercase.")
    if not data.obs_names.is_unique:
        raise ValueError(f"{label} cell identifiers must be unique.")
    if data.X is None or min(data.shape) == 0:
        raise ValueError(f"{label} must contain a nonempty expression matrix in X.")
    values = data.X.data if sparse.issparse(data.X) else data.X
    if not np.isfinite(values).all():
        raise ValueError(f"{label} expression matrix contains nonfinite values.")
    return data


def read_genes(path):
    with open(path) as handle:
        return list(dict.fromkeys(line.strip().lower() for line in handle if line.strip()))


def predict(args):
    predicted = f"{args.method}_predicted_expression"
    spatial = read_data(args.spatial, "Spatial data")
    reference = read_data(args.reference, "Reference data")
    if args.method == "spage" and reference.n_obs < 50:
        raise ValueError("SpaGE requires at least 50 reference cells.")
    if args.spatial_key not in spatial.obsm:
        raise ValueError(f"Spatial coordinates are required in obsm[{args.spatial_key!r}].")
    coordinates = np.asarray(spatial.obsm[args.spatial_key])
    if coordinates.ndim != 2 or coordinates.shape[1] < 2 or not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates must be finite with at least two columns.")
    # TISSUE 1.0.1 accepts a spatial key but does not forward it to Squidpy.
    spatial.obsm["spatial"] = coordinates
    if args.targets:
        targets = read_genes(args.targets)
    else:
        targets = list(reference.var_names[~reference.var_names.isin(spatial.var_names)])
    if not targets:
        raise ValueError("No target genes found. Supply a nonempty gene list or a reference with unmeasured genes.")
    missing = sorted(set(targets) - set(reference.var_names))
    if missing:
        raise ValueError("Target genes absent from reference: " + ", ".join(missing))

    # Withhold explicitly requested measured genes to avoid training leakage.
    shared = np.setdiff1d(np.intersect1d(spatial.var_names, reference.var_names), targets)
    spatial = spatial[:, shared].copy()
    calibration_genes = read_genes(args.conf_genes) if args.conf_genes else list(shared)
    if not calibration_genes:
        raise ValueError("No calibration genes found. Supply shared, non-target genes.")
    missing = sorted(set(calibration_genes) - set(shared))
    if missing:
        raise ValueError("Calibration genes must be shared, non-target genes: " + ", ".join(missing))
    if args.n_folds > len(calibration_genes):
        raise ValueError("Cross-validation folds cannot exceed the number of calibration genes.")
    smallest_training_set = len(shared) - int(np.ceil(len(calibration_genes) / args.n_folds))
    if args.method == "spage" and args.n_pv > min(smallest_training_set, spatial.n_obs, reference.n_obs):
        raise ValueError("Principal vectors must fit the smallest cross-validation training set and both cell counts.")
    if args.method == "knn":
        if args.knn_neighbors > reference.n_obs:
            raise ValueError("Imputation neighbors cannot exceed the reference cell count.")
        if smallest_training_set < 2:
            raise ValueError("Harmony-kNN requires at least two shared training genes in every fold.")
    if args.graph_method in ("fixed", "fixed_radius") and args.n_neighbors >= spatial.n_obs:
        raise ValueError("Spatial neighbors must be smaller than the spatial cell count.")

    if args.normalize_reference:
        values = reference.X.data if sparse.issparse(reference.X) else reference.X
        if (values < 0).any() or (np.asarray(reference.X.sum(axis=1)).ravel() <= 0).any():
            raise ValueError("Reference normalization requires nonnegative counts and positive cell totals.")
        tissue.main.preprocess_data(reference, standardize=False, normalize=True)

    # TISSUE's calibration residual calculation expects a dense spatial matrix.
    if sparse.issparse(spatial.X):
        spatial.X = spatial.X.toarray()
    if args.method == "spage":
        for data, label in ((spatial, "Spatial data"), (reference[:, shared], "Reference data")):
            matrix = data.X.toarray() if sparse.issparse(data.X) else data.X
            if (np.std(matrix, axis=0) == 0).any():
                raise ValueError(f"{label} contains constant shared genes; remove them before running SpaGE.")
    # Upstream Harmony-kNN stacks X with numpy and requires dense input.
    if args.method == "knn" and sparse.issparse(reference.X):
        reference.X = reference.X.toarray()

    np.random.seed(args.seed)
    if args.method == "tangram":
        import torch

        torch.manual_seed(args.seed)
    method_options = {"n_pv": args.n_pv} if args.method == "spage" else {}
    if args.method == "knn":
        method_options["n_neighbors"] = args.knn_neighbors
    tissue.main.predict_gene_expression(
        spatial, reference, targets,
        conf_genes=calibration_genes.copy() if args.conf_genes else None,
        method=args.method, n_folds=args.n_folds,
        random_seed=args.seed, **method_options,
    )
    predictions = spatial.obsm[predicted]
    if not np.isfinite(predictions.values).all():
        raise ValueError(f"{args.method} produced nonfinite predictions; check expression variability and method settings.")
    if args.gene_groups > predictions.shape[1] or args.cell_groups > spatial.n_obs:
        raise ValueError("Group counts cannot exceed the predicted gene count or spatial cell count.")
    if (args.gene_groups == 0 and predictions.shape[1] < 5) or (args.cell_groups == 0 and spatial.n_obs < 5):
        raise ValueError("Automatic grouping requires at least five genes or cells, respectively.")
    if args.weight == "exp_cos" and args.weight_n_pc > min(predictions.shape):
        raise ValueError("Weighting PCs cannot exceed either dimension of the prediction matrix (0 uses all genes).")

    tissue.main.build_spatial_graph(
        spatial, method=args.graph_method, n_neighbors=args.n_neighbors,
        radius=args.radius, set_diag=True,
    )
    tissue.main.conformalize_spatial_uncertainty(
        spatial, predicted, calib_genes=calibration_genes,
        grouping_method="kmeans_gene_cell", k=args.gene_groups or "auto",
        k2=args.cell_groups or "auto", weight=args.weight,
        weight_n_pc=args.weight_n_pc or None, add_one=True,
    )
    tissue.main.conformalize_prediction_interval(
        spatial, predicted, calib_genes=calibration_genes,
        alpha_level=args.alpha, symmetric=args.symmetric,
        compute_wasserstein=args.compute_wasserstein,
    )
    for suffix, path in (("", args.predictions), ("_lo", args.lower),
                         ("_hi", args.upper), ("_uncertainty", args.uncertainty),
                         ("_wasserstein", args.wasserstein)):
        if path:
            spatial.obsm[predicted + suffix].loc[:, targets].to_csv(
                path, sep="\t", index_label="cell_id", float_format="%.10g", na_rep="NA",
            )
    spatial.uns["tissue_galaxy"] = {
        "tissue_version": version("tissue-sc"), "seed": args.seed,
        "method": args.method,
        "n_folds": args.n_folds,
        "normalize_reference": args.normalize_reference,
        "symmetric": args.symmetric,
    }
    if args.method == "spage":
        spatial.uns["tissue_galaxy"]["n_pv"] = args.n_pv
    elif args.method == "knn":
        spatial.uns["tissue_galaxy"]["knn_neighbors"] = args.knn_neighbors
    spatial.write_h5ad(args.output, compression="gzip")
    print(f"Predicted {len(targets)} target genes in {spatial.n_obs} cells; calibrated with {len(calibration_genes)} shared genes.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=version("tissue-sc"))
    parser.add_argument("--spatial", required=True)
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--targets")
    parser.add_argument("--conf-genes")
    parser.add_argument("--normalize-reference", action="store_true")
    parser.add_argument("--method", choices=["spage", "knn", "tangram"], default="spage")
    parser.add_argument("--knn-neighbors", type=int, default=10)
    parser.add_argument("--n-folds", type=int, default=10)
    parser.add_argument("--n-pv", type=int, default=10)
    parser.add_argument("--seed", type=int, default=444)
    parser.add_argument("--graph-method", choices=["fixed_radius", "fixed", "delaunay", "delaunay_radius", "radius"], default="fixed_radius")
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--radius", type=float)
    parser.add_argument("--gene-groups", type=int, default=0)
    parser.add_argument("--cell-groups", type=int, default=0)
    parser.add_argument("--weight", choices=["exp_cos", "spatial_connectivities"], default="exp_cos")
    parser.add_argument("--weight-n-pc", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.33)
    parser.add_argument("--symmetric", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--compute-wasserstein", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--lower", required=True)
    parser.add_argument("--upper", required=True)
    parser.add_argument("--uncertainty", required=True)
    parser.add_argument("--wasserstein")
    args = parser.parse_args()
    if args.compute_wasserstein != bool(args.wasserstein):
        parser.error("--compute-wasserstein and --wasserstein must be supplied together")
    try:
        predict(args)
    except ValueError as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
