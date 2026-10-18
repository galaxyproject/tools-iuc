#!/usr/bin/env python
"""Run the TISSUE Python API and export its results for Galaxy."""

import argparse
from importlib.metadata import version

import anndata
import numpy as np
import tissue.main
from scipy import sparse

__version__ = "1.0.1"
PREDICTED = "spage_predicted_expression"


def read_data(path, label):
    data = anndata.read_h5ad(path)
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


def predict(args):
    spatial = read_data(args.spatial, "Spatial data")
    reference = read_data(args.reference, "Reference data")
    if reference.n_obs < 50:
        raise ValueError("SpaGE requires at least 50 reference cells.")
    if "spatial" not in spatial.obsm:
        raise ValueError("Spatial coordinates are required in obsm['spatial'].")
    coordinates = np.asarray(spatial.obsm["spatial"])
    if coordinates.ndim != 2 or coordinates.shape[1] < 2 or not np.isfinite(coordinates).all():
        raise ValueError("Spatial coordinates must be finite with at least two columns.")
    if args.targets:
        with open(args.targets) as handle:
            targets = list(dict.fromkeys(line.strip().lower() for line in handle if line.strip()))
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
    if args.n_folds < 2 or args.n_folds > len(shared):
        raise ValueError("Cross-validation folds must be between 2 and the number of shared, non-target genes.")
    smallest_training_set = len(shared) - int(np.ceil(len(shared) / args.n_folds))
    if not 1 <= args.n_pv <= min(smallest_training_set, spatial.n_obs, reference.n_obs):
        raise ValueError("Principal vectors must fit the smallest cross-validation training set and both cell counts.")
    if args.graph_method in ("fixed", "fixed_radius") and not 1 <= args.n_neighbors < spatial.n_obs:
        raise ValueError("Spatial neighbors must be positive and smaller than the spatial cell count.")
    if args.radius is not None and args.radius <= 0:
        raise ValueError("Spatial radius must be positive.")

    if args.normalize_reference:
        values = reference.X.data if sparse.issparse(reference.X) else reference.X
        if (values < 0).any() or (np.asarray(reference.X.sum(axis=1)).ravel() <= 0).any():
            raise ValueError("Reference normalization requires nonnegative counts and positive cell totals.")
        tissue.main.preprocess_data(reference, standardize=False, normalize=True)

    # TISSUE's calibration residual calculation expects a dense spatial matrix.
    if sparse.issparse(spatial.X):
        spatial.X = spatial.X.toarray()
    for data, label in ((spatial, "Spatial data"), (reference[:, shared], "Reference data")):
        matrix = data.X.toarray() if sparse.issparse(data.X) else data.X
        if (np.std(matrix, axis=0) == 0).any():
            raise ValueError(f"{label} contains constant shared genes; remove them before running SpaGE.")

    np.random.seed(args.seed)
    tissue.main.predict_gene_expression(
        spatial, reference, targets, method="spage", n_folds=args.n_folds,
        random_seed=args.seed, n_pv=args.n_pv,
    )
    predictions = spatial.obsm[PREDICTED]
    if not np.isfinite(predictions.values).all():
        raise ValueError("SpaGE produced nonfinite predictions; check expression variability and principal vectors.")
    if args.gene_groups < 0 or args.cell_groups < 0:
        raise ValueError("Group counts must be nonnegative (0 selects automatic grouping).")
    if args.gene_groups > predictions.shape[1] or args.cell_groups > spatial.n_obs:
        raise ValueError("Group counts cannot exceed the predicted gene count or spatial cell count.")
    if (args.gene_groups == 0 and predictions.shape[1] < 5) or (args.cell_groups == 0 and spatial.n_obs < 5):
        raise ValueError("Automatic grouping requires at least five genes or cells, respectively.")
    if args.weight == "exp_cos" and not 0 <= args.weight_n_pc <= min(predictions.shape):
        raise ValueError("Weighting PCs cannot exceed either dimension of the prediction matrix (0 uses all genes).")
    if not 0 < args.alpha < 1:
        raise ValueError("Alpha must be strictly between 0 and 1.")

    tissue.main.build_spatial_graph(
        spatial, method=args.graph_method, n_neighbors=args.n_neighbors,
        radius=args.radius, set_diag=True,
    )
    tissue.main.conformalize_spatial_uncertainty(
        spatial, PREDICTED, calib_genes=spatial.var_names,
        grouping_method="kmeans_gene_cell", k=args.gene_groups or "auto",
        k2=args.cell_groups or "auto", weight=args.weight,
        weight_n_pc=args.weight_n_pc or None, add_one=True,
    )
    tissue.main.conformalize_prediction_interval(
        spatial, PREDICTED, calib_genes=spatial.var_names,
        alpha_level=args.alpha, symmetric=not args.asymmetric,
        compute_wasserstein=args.compute_wasserstein,
    )
    for suffix, path in (("", args.predictions), ("_lo", args.lower),
                         ("_hi", args.upper), ("_uncertainty", args.uncertainty),
                         ("_wasserstein", args.wasserstein)):
        if path:
            spatial.obsm[PREDICTED + suffix].loc[:, targets].to_csv(
                path, sep="\t", index_label="cell_id", float_format="%.10g", na_rep="NA",
            )
    spatial.uns["tissue_galaxy"] = {
        "tissue_version": version("tissue-sc"), "seed": args.seed,
        "n_folds": args.n_folds, "n_pv": args.n_pv,
        "normalize_reference": args.normalize_reference,
        "symmetric": not args.asymmetric,
    }
    spatial.write_h5ad(args.output, compression="gzip")
    print(f"Predicted {len(targets)} target genes in {spatial.n_obs} cells; calibrated with {len(shared)} shared genes.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=version("tissue-sc"))
    parser.add_argument("--spatial", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--targets")
    parser.add_argument("--normalize-reference", action="store_true")
    parser.add_argument("--n-folds", type=int, default=10)
    parser.add_argument("--n-pv", type=int, default=10)
    parser.add_argument("--seed", type=int, default=444)
    parser.add_argument("--graph-method", choices=["fixed_radius", "fixed", "delaunay", "delaunay_radius", "radius"], default="fixed_radius")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--radius", type=float)
    parser.add_argument("--gene-groups", type=int, default=4)
    parser.add_argument("--cell-groups", type=int, default=2)
    parser.add_argument("--weight", choices=["exp_cos", "spatial_connectivities"], default="exp_cos")
    parser.add_argument("--weight-n-pc", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.33)
    parser.add_argument("--asymmetric", action="store_true")
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
