The different methods from Scanpy have been grouped by themes:

1. Filter in `filter.xml`
  - Filter cell outliers based on counts and numbers of genes expressed, using `pp.filter_cells`
  - Filter genes based on number of cells or counts, using `pp.filter_genes`
  - Extract highly variable genes, using `pp.filter_genes_dispersion`
  - `tl.highly_variable_genes` (need to be added)
  - Subsample to a fraction of the number of observations, using `pp.subsample`
  - `queries.gene_coordinates` (need to be added)
  - `queries.mitochondrial_genes` (need to be added)

2. Normalize in `normalize.xml`
  - Normalize total counts per cell, using `pp.normalize_per_cell`
  - Normalization and filtering as of Zheng et al. (2017), using `pp.recipe_zheng17`
  - Normalization and filtering as of Weinreb et al (2017), using `pp.recipe_weinreb17`
  - Normalization and filtering as of Seurat et al (2015), using `pp.recipe_seurat`
  - Logarithmize the data matrix, using `pp.log1p`
  - Scale data to unit variance and zero mean, using `pp.scale`
  - Square root the data matrix, using `pp.sqrt`
  - Downsample counts, using `pp.downsample_counts`

3. Remove confounder in `remove_confounders.xml`
  - Regress out unwanted sources of variation, using `pp.regress_out`
  - `pp.mnn_correct` (need to be added)
  - `pp.mnn_correct` (need to be added)
  - `pp.magic` (need to be added)
  - `tl.sim` (need to be added)
  - `pp.calculate_qc_metrics` (need to be added)
  - Score a set of genes, using `tl.score_genes`
  - Score cell cycle genes, using `tl.score_genes_cell_cycle`
  - `tl.cyclone` (need to be added)
  - `tl.andbag` (need to be added)

4. Cluster and reduce dimension in `cluster_reduce_dimension.xml`
  - `tl.leiden` (need to be added)
  - Cluster cells into subgroups, using `tl.louvain`
  - Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using `pp.pca`
  - Computes PCA (principal component analysis) coordinates, loadings and variance decomposition, using `tl.pca`
  - Diffusion Maps, using `tl.diffmap`
  - t-distributed stochastic neighborhood embedding (tSNE), using `tl.tsne`
  - Embed the neighborhood graph using UMAP, using `tl.umap`
  - `tl.phate` (need to be added)
  - Compute a neighborhood graph of observations, using `pp.neighbors`
  - Rank genes for characterizing groups, using `tl.rank_genes_groups`

4. Inspect
  - `tl.paga_compare_paths` (need to be added)
  - `tl.paga_degrees` (need to be added)
  - `tl.paga_expression_entropies` (need to be added)
  - Generate cellular maps of differentiation manifolds with complex topologies, using `tl.paga`
  - Infer progression of cells through geodesic distance along the graph, using `tl.dpt`

5. Plot
  1. Generic
    - Scatter plot along observations or variables axes, using `pl.scatter`
    - Heatmap of the expression values of set of genes, using `pl.heatmap`
    - Makes a dot plot of the expression values, using `pl.dotplot`
    - Violin plot, using `pl.violin`
    - `pl.stacked_violin` (need to be added)
    - Heatmap of the mean expression values per cluster, using `pl.matrixplot`
    - Hierarchically-clustered heatmap, using `pl.clustermap`
    - `pl.ranking` 

  2. Preprocessing
    - Plot the fraction of counts assigned to each gene over all cells, using `pl.highest_expr_genes`
    - Plot dispersions versus means for genes, using `pl.filter_genes_dispersion`
    - `pl.highly_variable_genes` (need to be added)
    - `pl.calculate_qc_metrics` (need to be added)
  
  3. PCA
    - Scatter plot in PCA coordinates, using `pl.pca`
    - Rank genes according to contributions to PCs, using `pl.pca_loadings`
    - Scatter plot in PCA coordinates, using `pl.pca_variance_ratio`
    - Plot PCA results, using `pl.pca_overview`
  
  4. Embeddings
    - Scatter plot in tSNE basis, using `pl.tsne`
    - Scatter plot in UMAP basis, using `pl.umap`
    - Scatter plot in Diffusion Map basis, using `pl.diffmap`
    - `pl.draw_graph` (need to be added)

  5. Branching trajectories and pseudotime, clustering
    - Plot groups and pseudotime, using `pl.dpt_groups_pseudotime`
    - Heatmap of pseudotime series, using `pl.dpt_timeseries`
    - Plot the abstracted graph through thresholding low-connectivity edges, using `pl.paga`
    - `pl.paga_compare` (need to be added)
    - `pl.paga_path` (need to be added)

  6. Marker genes: 
    - Plot ranking of genes using dotplot plot, using `pl.rank_gene_groups`
    - `pl.rank_genes_groups_dotplot` (need to be added)
    - `pl.rank_genes_groups_heatmap` (need to be added)
    - `pl.rank_genes_groups_matrixplot` (need to be added)
    - `pl.rank_genes_groups_stacked_violin` (need to be added)
    - `pl.rank_genes_groups_violin` (need to be added)

  7. Misc
    - `pl.phate` (need to be added)
    - `pl.matrix` (need to be added)
    - `pl.paga_adjacency` (need to be added)
    - `pl.timeseries` (need to be added)
    - `pl.timeseries_as_heatmap` (need to be added)
    - `pl.timeseries_subplot` (need to be added)
    
  