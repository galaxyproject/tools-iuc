Scanpy
======

1. Inspect & Manipulate (`inspect.xml`)

    Methods | Description
    --- | ---
    `pp.calculate_qc_metrics` | Calculate quality control metrics
    `pp.neighbors` | Compute a neighborhood graph of observations
    `tl.score_genes` | Score a set of genes
    `tl.score_genes_cell_cycle` | Score cell cycle gene
    `tl.rank_genes_groups` | Rank genes for characterizing groups
    `tl.marker_gene_overlap` | Calculate an overlap score between data-deriven marker genes and provided markers (**not working for now**)
    `pp.log1p` | Logarithmize the data matrix.
    `pp.scale` | Scale data to unit variance and zero mean
    `pp.sqrt` | Square root the data matrix

2. Filter (`filter.xml`)

    Methods | Description
    --- | ---
    `pp.filter_cells` | Filter cell outliers based on counts and numbers of genes expressed.
    `pp.filter_genes` | Filter genes based on number of cells or counts.
    `tl.filter_rank_genes_groups` | Filters out genes based on fold change and fraction of genes expressing the gene within and outside the groupby categories (**to fix**)
    `pp.highly_variable_genes` | Extract highly variable genes
    `pp.subsample` | Subsample to a fraction of the number of observations
    `pp.downsample_counts` | Downsample counts so that each cell has no more than target_counts

3. Normalize (`normalize.xml`)

    Methods | Description
    --- | ---
    `pp.normalize_total` | Normalize counts per cell
    `pp.recipe_zheng17` | Normalization and filtering as of [Zheng17]
    `pp.recipe_weinreb17` | Normalization and filtering as of [Weinreb17]
    `pp.recipe_seurat` | Normalization and filtering as of Seurat [Satija15]

4. Remove confounders (`remove_confounder.xml`)

    Methods | Description
    --- | ---
   `pp.regress_out` | Regress out unwanted sources of variation
   `pp.mnn_correct` | Correct batch effects by matching mutual nearest neighbors
   `pp.combat` | ComBat function for batch effect correction

5. Clustering, embedding and trajectory inference (`cluster_reduce_dimension.xml`)

    Methods | Description
    --- | ---
    `tl.louvain` | Cluster cells into subgroups
    `tl.leiden` | Cluster cells into subgroups
    `tl.pca` | Principal component analysis
    `pp.pca` | Principal component analysis (appears to be the same func...)
    `tl.diffmap` | Diffusion Maps
    `tl.tsne` | t-SNE
    `tl.umap` | Embed the neighborhood graph using UMAP
    `tl.draw_graph` | Force-directed graph drawing
    `tl.dpt` | Infer progression of cells through geodesic distance along the graph
    `tl.paga` | Mapping out the coarse-grained connectivity structures of complex manifolds

6. Plot (`plot.xml`)

    1. Generic

        Methods | Description
        --- | ---
        `pl.scatter` | Scatter plot along observations or variables axes
        `pl.heatmap` | Heatmap of the expression values of set of genes
        `pl.dotplot` | Makes a dot plot of the expression values
        `pl.violin` | Violin plot
        `pl.stacked_violin` | Stacked violin plots
        `pl.matrixplot` | Heatmap of the mean expression values per cluster
        `pl.clustermap` | Hierarchically-clustered heatmap
    
    2. Preprocessing

        Methods | Description
        --- | ---
        `pl.highest_expr_genes` | Plot the fraction of counts assigned to each gene over all cells
        `pl.highly_variable_genes` | Plot dispersions versus means for genes

    3. PCA

        Methods | Description
        --- | ---
        `pl.pca` | Scatter plot in PCA coordinates
        `pl.pca_loadings` | Rank genes according to contributions to PCs
        `pl.pca_variance_ratio` | Scatter plot in PCA coordinates
        `pl.pca_overview` | Plot PCA results

    4. Embeddings

        Methods | Description
        --- | ---
        `pl.tsne` | Scatter plot in tSNE basis
        `pl.umap` | Scatter plot in UMAP basis
        `pl.diffmap` | Scatter plot in Diffusion Map basis
        `pl.draw_graph` | Scatter plot in graph-drawing basis

    5. Branching trajectories and pseudotime, clustering

        Methods | Description
        --- | ---
        `pl.dpt_groups_pseudotime` | Plot groups and pseudotime
        `pl.dpt_timeseries` | Heatmap of pseudotime series
        `pl.paga` | Plot the abstracted graph through thresholding low-connectivity edges
        `pl.paga_compare` | Scatter and PAGA graph side-by-side
        `pl.paga_path` | Gene expression and annotation changes along paths

    6. Marker genes

        Methods | Description
        --- | ---
        `pl.rank_genes_groups` | Plot ranking of genes using dotplot plot
        `pl.rank_genes_groups_violin` | Plot ranking of genes for all tested comparisons
