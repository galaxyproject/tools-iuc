Scanpy
======

## Classification of methods into steps

Steps:
1. Filtering

    Methods | Description
    --- | ---
    `pp.filter_cells` | Filter cell outliers based on counts and numbers of genes expressed.
    `pp.filter_genes` | Filter genes based on number of cells or counts.
    `pp.filter_genes_dispersion` | Extract highly variable genes
    `pp.highly_variable_genes` | Extract highly variable genes
    `pp.subsample` | Subsample to a fraction of the number of observations
    `queries.gene_coordinates` | (Could not find...)
    `queries.mitochondrial_genes` | Retrieves Mitochondrial gene symbols for specific organism through BioMart for filtering

2. Quality Plots

   These are in-between stages used to measure the effectiveness of a Filtering/Normalisation/Conf.Removal stage either after processing or prior to.

    Methods | Description | Notes
    --- | --- | ---
    `pp.calculate_qc_metrics` | Calculate quality control metrics
    `pl.violin` | violin plot of features, lib. size, or subsets of. 
    `pl.stacked_violin` | Same as above but for multiple series of features or cells

3. Normalization

    Methods | Description
    --- | ---
    `pp.normalize_per_cell` | Normalize total counts per cell
    `pp.recipe_zheng17` | Normalization and filtering as of [Zheng17]
    `pp.recipe_weinreb17` | Normalization and filtering as of [Weinreb17]
    `pp.recipe_seurat` | Normalization and filtering as of Seurat [Satija15]
    `pp.log1p` | Logarithmize the data matrix.
    `pp.scale` | Scale data to unit variance and zero mean
    `pp.sqrt` | 
    `pp.downsample_counts` | Downsample counts so that each cell has no more than target_counts

4. Conf. removal

    Methods | Description
    --- | ---
   `pp.regress_out` | Regress out unwanted sources of variation
   `pp.mnn_correct` | Correct batch effects by matching mutual nearest neighbors
   `pp.dca` | Deep count autoencoder to denoise the data
   `pp.magic` | Markov Affinity-based Graph Imputation of Cells (MAGIC) API to denoise
   `tl.sim` | Simulate dynamic gene expression data [Wittman09]
   `pp.calculate_qc_metrics` | Calculate quality control metrics
   `tl.score_genes` | Score a set of genes
   `tl.score_genes_cell_cycle` | Score cell cycle genes
   `tl.cyclone` | Assigns scores and predicted class to observations based on cell-cycle genes [Scialdone15]
   `tl.sandbag` | Calculates pairs of genes serving as markers for each cell-cycle phase [Scialdone15]

5. Clustering and Heatmaps

    Methods | Description
    --- | ---
    `tl.leiden` | Cluster cells into subgroups [Traag18] [Levine15]
    `tl.louvain` | Cluster cells into subgroups [Blondel08] [Levine15] [Traag17]
    `tl.pca` | Principal component analysis
    `pp.pca` | Principal component analysis (appears to be the same func...)
    `tl.diffmap` | Diffusion Maps
    `tl.tsne` | t-SNE
    `tl.umap` | Embed the neighborhood graph using UMAP
    `tl.phate` | PHATE
    `pp.neighbors` | Compute a neighborhood graph of observations
    `tl.rank_genes_groups` | Rank genes for characterizing groups
    `pl.rank_genes_groups` | 
    `pl.rank_genes_groups_dotplot` | 
    `pl.rank_genes_groups_heatmap` | 
    `pl.rank_genes_groups_matrixplot` | 
    `pl.rank_genes_groups_stacked_violin` | 
    `pl.rank_genes_groups_violin` | 
    `pl.matrix_plot` | 
    `pl.heatmap` | 
    `pl.highest_expr_genes` | 
    `pl.diffmap` | 
    
6. Cluster Inspection and plotting

    Methods that draw out the clusters computed in the previous stage, not heatmap or pseudotime related.

    Methods | Description 
    --- | --- 
    `pl.clustermap` |
    `pl.phate` | 
    `pl.dotplot` | 
    `pl.draw_graph` | (really general purpose, would not implement directly)
    `pl.filter_genes_dispersion` | (depreciated for 'highly_variable_genes')
    `pl.matrix` | (could not find in API)
    `pl.pca` | 
    `pl.pca_loadings` | 
    `pl.pca_overview` | 
    `pl.pca_variance_ratio` | 
    `pl.ranking` | (not sure what this does...)
    `pl.scatter` | ([very general purpose](https://icb-scanpy.readthedocs-hosted.com/en/latest/api/scanpy.api.pl.scatter.html), would not implement directly)
    `pl.set_rcParams_defaults` | 
    `pl.set_rcParams_scanpy` | 
    `pl.sim` | 
    `pl.tsne` | 
    `pl.umap` | 

7. Branch/Between-Cluster Inspection

    Pseudotime analysis, relies on initial clustering.

    Methods | Description
    --- | ---
    `tl.dpt` | Infer progression of cells through geodesic distance along the graph [Haghverdi16] [Wolf17i]
    `pl.dpt_groups_pseudotime` | 
    `pl.dpt_timeseries` | 
    `tl.paga_compare_paths` | 
    `tl.paga_degrees` | 
    `tl.paga_expression_entropies` | 
    `tl.paga` | Generate cellular maps of differentiation manifolds with complex topologies [Wolf17i]
    `pl.paga` | 
    `pl.paga_adjacency` | 
    `pl.paga_compare` | 
    `pl.paga_path` | 
    `pl.timeseries` | 
    `pl.timeseries_as_heatmap` | 
    `pl.timeseries_subplot` | 


Methods to sort | Description
--- | --- 
`tl.ROC_AUC_analysis` | (could not find in API)
`tl.correlation_matrix` | (could not find in API)
`rtools.mnn_concatenate` | (could not find in API)
`utils.compute_association_matrix_of_groups` | (could not find in API) 
`utils.cross_entropy_neighbors_in_rep` | (could not find in API)
`utils.merge_groups` | (could not find in API)
`utils.plot_category_association` | (could not find in API)
`utils.select_groups` | (could not find in API)