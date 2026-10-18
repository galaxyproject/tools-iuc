# TEST FILE FOR BASIC TISSUE FUNCTIONALITIES


# import packages

import tissue.main, tissue.downstream

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import os

import warnings
warnings.filterwarnings("ignore")

#################################################################################################################
print ("Testing TISSUE data loading...")
try:
    adata, RNAseq_adata = tissue.main.load_paired_datasets("tests/data/Spatial_count.txt",
                                                           "tests/data/Locations.txt",
                                                           "tests/data/scRNA_count.txt")
except:
    raise Exception ("Failed data loading from tests/data/ with tissue.main.load_paired_datasets()")

#################################################################################################################
print ("Testing TISSUE preprocessing...")
adata.var_names = [x.lower() for x in adata.var_names]
RNAseq_adata.var_names = [x.lower() for x in RNAseq_adata.var_names]
try:
    tissue.main.preprocess_data(RNAseq_adata, standardize=False, normalize=True)
except:
    raise Exception ("Failed TISSUE preprocessing. Make sure all dependencies are installed.")
gene_names = np.intersect1d(adata.var_names, RNAseq_adata.var_names)
adata = adata[:, gene_names].copy()
target_gene = "plp1"
target_expn = adata[:, target_gene].X.copy()
adata = adata[:, [gene for gene in gene_names if gene != target_gene]].copy()

#################################################################################################################
print("Testing TISSUE spatial gene expression prediction...")
try:
    tissue.main.predict_gene_expression (adata, RNAseq_adata, [target_gene],
                                         method="spage", n_folds=3, n_pv=10)
except:
    raise Exception("TISSUE prediction failed for SpaGE at tissue.main.predict_gene_expression()")

#################################################################################################################
print("Testing TISSUE calibration...")
try:
    tissue.main.build_spatial_graph(adata, method="fixed_radius", n_neighbors=15)
except:
    raise Exception ("Failed TISSUE spatial graph building at tissue.main.build_spatial_graph()")
try:
    tissue.main.conformalize_spatial_uncertainty(adata, "spage_predicted_expression", calib_genes=adata.var_names,
                                                 grouping_method="kmeans_gene_cell", k=4, k2=2)
except:
    raise Exception ("Failed TISSUE cell-centric variability and calibration scores processing at tissue.main.conformalize_spatial_uncertainty()")
try:
    tissue.main.conformalize_prediction_interval (adata, "spage_predicted_expression", calib_genes=adata.var_names,
                                                  alpha_level=0.23, compute_wasserstein=True)
except:
    raise Exception ("Failed TISSUE prediction interval calibration at tissue.main.conformalize_prediction_interval()")

#################################################################################################################
print ("Testing TISSUE multiple imputation t-test...")
adata.obs['condition'] = ['A' if i < round(adata.shape[0]/2) else 'B' for i in range(adata.shape[0])]
try:
    tissue.downstream.multiple_imputation_testing(adata, "spage_predicted_expression",
                                                  calib_genes=adata.var_names,
                                                  condition='condition',
                                                  group1 = "A", # use None to compute for all conditions, condition vs all
                                                  group2 = "B", # use None to compute for group1 vs all
                                                  n_imputations=2)
except:
    raise Exception ("Failed TISSUE MI t-test at tissue.downstream.multiple_imputation_testing()")

#################################################################################################################
print("Testing TISSUE cell filtering")
X_uncertainty = adata.obsm["spage_predicted_expression_hi"].values - adata.obsm["spage_predicted_expression_lo"].values
try:
    keep_idxs = tissue.downstream.detect_uncertain_cells (X_uncertainty,
                                                          proportion="otsu",
                                                          stratification=adata.obs['condition'].values)
except:
    raise Exception ("Failed TISSUE cell filtering at tissue.downstream.detect_uncertain_cells()")
try:
    keep_idxs = tissue.downstream.filtered_PCA (adata, # anndata object
                                                "spage", # prediction method
                                                proportion="otsu",
                                                stratification=adata.obs['condition'].values,
                                                return_keep_idxs=True)
except:
    raise Exception ("Failed TISSUE-filtered PCA at tissue.downstream.filtered_PCA()")

print("TISSUE tests passed!")