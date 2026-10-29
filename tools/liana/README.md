LIANA+
======

Galaxy tools for [LIANA+](https://liana.readthedocs.io/en/stable/), a framework for cell-cell communication analysis. These wrappers target LIANA **2.0.0**.

## Tools

### 1. LIANA Methods (`liana_methods.xml`)

The combined methods tool contains single-cell ligand-receptor inference and spatial methods.

Method | Description
--- | ---
`cellchat` | CellChat ligand-receptor method
`cellphonedb` | CellPhoneDB ligand-receptor method
`connectome` | Connectome ligand-receptor method
`logfc` | Log fold-change ligand-receptor method
`natmi` | NATMI ligand-receptor method
`singlecellsignalr` | SingleCellSignalR ligand-receptor method
`geometric_mean` | Geometric-mean ligand-receptor method
`scseqcomm` | Magnitude-only ligand-receptor scoring from expression CDFs
`rank_aggregate` | Aggregate rankings from multiple methods
`bivariate` | Local/global bivariate spatial statistics for AnnData or paired MuData modalities; one or both global metrics can be selected
`cross_pcf` | Distance-resolved cross pair-correlation for unordered cell-type pairs
`lric` | Expression-weighted Ligand-Receptor Interaction Correlation
`compute_global_specificity` | Permutation-tested group specificity for global interactions
`inflow` | Trivariate spatial source-group, ligand, and receptor scoring

`cross_pcf` and `lric` return both an updated AnnData file and a long-form spatial-curves table suitable for `lric_lineplot`.

Single-cell ligand-receptor methods, including `rank_aggregate`, compare groups defined by an observation column. They do not use spatial coordinates or a spatial-neighbor graph to restrict the inferred interactions. For spatial analyses, select `bivariate`, `cross_pcf`, `lric`, or `inflow` according to the intended statistic.

Choose an expression source that satisfies the selected method's requirements. Single-cell ligand-receptor inference requires nonnegative, log-normalized expression; a scaled or centered matrix is unsuitable. Enabling **Use adata.raw** requires an initialized `.raw` matrix. It does not fall back automatically to `.X`. To use `.X` or a named layer, disable that option. Cell-identity or predicted-label columns remain separate annotations: the selected **Group by** column determines the groups analyzed.

### 2. MISTy (`misty.xml`)

Method | Description
--- | ---
`genericMistyData` | Generic MISTy multi-view modelling
`lrMistyData` | Ligand-receptor MISTy modelling

MISTy accepts AnnData-derived views or an existing MuData object, depending on the selected input mode. The annotated output contains fitted model results; plot the target-metric and contribution tables with the matching MISTy plotting methods.

### 3. Multi-view Utilities (`multi.xml`)

Wrappers for converting LIANA results to views/tensors and for NMF-based multi-view analysis. The Galaxy forms expose typed controls for `adata_to_views` pseudobulk and filtering options, `lrs_to_views` batch-aware variance filtering and missing-value filling, and stable NMF initialization/convergence controls.

Method | Description
--- | ---
`to_tensor_c2c` | Convert LIANA results to a cell2cell-compatible tensor
`lrs_to_views` | Convert sample-level interaction results into MuData views
`adata_to_views` | Build cell-type pseudobulk MuData views
`df_to_lr` | Combine differential-expression statistics with an interaction resource
`lrdata_to_mudata` | Split LIANA feature data into paired MuData modalities
`filter_view_markers` | Flag or retain view-specific marker features
`nmf` | Non-negative matrix factorization with inclusive rank search bounds

### 4. Plotting (`plot.xml`)

Method | Description
--- | ---
`dotplot` | Dot plot of ligand-receptor interactions
`dotplot_by_sample` | Dot plot grouped by sample
`tileplot` | Tile plot of interactions
`connectivity` | Spatial connectivity visualization
`target_metrics` | MISTy target metrics
`contributions` | MISTy view contributions
`interactions` | MISTy interaction importances
`annulus_plot` | Inspect the spatial annulus geometry used by cross-PCF/LRIC
`lric_lineplot` | Plot distance-resolved g(r) curves from the methods output table
`circle_plot` | Aggregate directed interactions as a circular network
`feature_by_group` | Plot a spatial feature in group-specific panels

### 5. Utilities (`utils.xml`)

Includes spatial-neighbor construction, AnnData extraction helpers, factor/loadings extraction, `expand_coordinates`, MuData-to-AnnData conversion, spatial pair-proximity summaries, spatial interpolation between reference and target AnnData objects, `query_bandwidth` diagnostics with plot and table outputs, and `get_lric_auc` summaries of saved LRIC/cross-PCF curves.

The curve-summary utility reads the full spatial result from an AnnData uns key and writes interaction identifiers, the span-normalized AUC score, and peak radius. Its default log2 transformation floors g at 0.05: positive scores indicate co-enrichment, negative scores indicate depletion, and zero is the null. Distance limits are exclusive; curves with insufficient finite radius bins are omitted. Input expression and annotations are not changed.

### 6. Resources (`resource.xml`)

Provides built-in or Data Manager-cached ligand-receptor resources, Metalinks retrieval/filtering, HCOP ortholog retrieval through the current target-organism API, and resource translation.

## API coverage and Galaxy-specific constraints

The wrappers call LIANA 2.0.0 for the exposed computations. Python-native outputs are converted only after computation into Galaxy-compatible AnnData, MuData, tabular, or image outputs. Applicable LIANA methods receive `n_jobs` from Galaxy's allocated `GALAXY_SLOTS` value rather than from a user-controlled field. The public [API reference](https://liana.readthedocs.io/en/stable/api.html) describes the underlying methods; the options exposed by each Galaxy form are documented in its help panel.

The following Python-only extension points are intentionally not exposed:

- arbitrary callables such as `lrs_to_views(inverse_fun=...)`;
- callable or dictionary pseudobulk aggregation modes;
- unrestricted `**kwargs` or user-entered Python/JSON dictionaries.

Instead, stable options are represented with typed Galaxy parameters and conditionals. This keeps jobs reproducible, validates values before execution, and avoids evaluating user-supplied Python code.

The deterministic local fixtures combine numerical and structural assertions. Spatial reference scores use LIANA 2.0's random-labelling null and retain the complete long-format result, including `g_pcf` and `g_expr` for grouped LRIC. Other numerical assertions are retained where the expected behavior is unchanged.

## References

- Dimitrov, D., Türei, D., Garrido-Rodriguez, M. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0
- Tanevski, J., Flores, R.O.R., Gabor, A. et al. Explainable multiview framework for dissecting spatial relationships from highly multiplexed data. Genome Biol 23, 97 (2022).
- [LIANA+ Documentation](https://liana.readthedocs.io/en/stable/)
- [LIANA+ GitHub](https://github.com/scverse/liana)

## Data Manager integration

Galaxy administrators install and run `data_manager_liana_resources`. It registers ligand-receptor, HCOP, and Metalinks files in the shared `liana_resources` Tool Data Table. Regular users do not run the Data Manager: the LIANA tool forms expose those server-side entries through `from_data_table`.

For `translate_resource`, the ligand-receptor source and the ortholog mapping source each support a cached mode. The ortholog selector also supports `builtin`, which invokes `li.rs.get_hcop_orthologs` for the selected target organism, and `history`, which accepts a user-supplied table.

### Tool Data Table and tests

The installation templates are `tool_data_table_conf.xml.sample` and `tool-data/liana_resources.loc.sample`. The shared table has four columns: `value`, `name`, `description`, and `path`. Cached tool inputs use the registered path, rather than reconstructing a server-specific directory. For Planemo tests, `tool_data_table_conf.xml.test` points to `test-data/liana_resources.loc`, whose entries resolve local ligand-receptor, HCOP, and Metalinks fixtures through `${__HERE__}`.

Repeated spatial-annulus, cell-type, HCOP-download, MuData, separator, observation-key, view-key, and result-key parameters are defined in `macros.xml` so method categories share identical defaults and validation.

## LIANA 2.0 compatibility

The Galaxy method selectors are retained for existing workflows. Preprocessing calls now use `li.pp`, multisample calls use `li.ms`, and differential-expression interaction conversion uses `li.mt.df_to_lr`. These functions were relocated rather than removed. Tensor conversion uses cell2cell, expression-view conversion uses decoupler, and NMF rank selection uses kneed; these dependencies are supplied by the LIANA Bioconda package. Tensor export writes a numeric NumPy array without pickle and a separate table mapping each tensor axis index to its biological label.

Cross-PCF emits each unordered cell-type pair once. Grouped LRIC remains directed. Both retain every column returned by LIANA in the H5AD and the tabular output. Existing annulus widths must be a positive integer multiple of the radius step. Angular sampling is no longer used by the analytical null. A nonzero legacy minimum-expressing-cell count must be replaced explicitly with the appropriate expression-proportion filter; the wrapper does not guess that scientific choice. Proportion filters default to zero to retain the previous wrapper's unfiltered default.

Nullable-string serialization is explicitly enabled for all annotated-data outputs. The nullable PBMC regression fixture checks observation IDs, annotation missingness and the LIANA result table. PBMC expression tests use the nonnegative, log-normalized raw matrix rather than the scaled matrix in X. No expression values are clipped to satisfy the new validation.
