LIANA+
======

Galaxy tools for [LIANA+](https://liana-py.readthedocs.io/), a framework for cell-cell communication analysis. This wrapper targets LIANA **1.8.1**.

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
`rank_aggregate` | Aggregate rankings from multiple methods
`bivariate` | Local/global bivariate spatial statistics
`cross_pcf` | Distance-resolved cross pair-correlation for directed cell-type pairs
`lric` | Expression-weighted Ligand-Receptor Interaction Correlation

`cross_pcf` and `lric` return both an updated AnnData file and a long-form spatial-curves table suitable for `lric_lineplot`.

### 2. MISTy (`misty.xml`)

Method | Description
--- | ---
`genericMistyData` | Generic MISTy multi-view modelling
`lrMistyData` | Ligand-receptor MISTy modelling

LIANA 1.8.1 includes improved preservation of MuData metadata (`uns`, `obsm`, `varm`, `obsp`, and `varp`) during MISTyData round-trips.

### 3. Multi-view Utilities (`multi.xml`)

Wrappers for converting LIANA results to views/tensors and for NMF-based multi-view analysis.

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

### 5. Utilities (`utils.xml`)

Includes spatial-neighbor construction, AnnData extraction helpers, factor/loadings extraction, and LIANA 1.8.1's `expand_coordinates` for laying out multiple spatial samples on a non-overlapping grid.

### 6. Resources (`resource.xml`)

Provides built-in or Data Manager-cached ligand-receptor resources, Metalinks retrieval/filtering, HCOP ortholog retrieval through the current target-organism API, and resource translation.

## References

- Dimitrov, D., Türei, D., Garrido-Rodriguez, M. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0
- Tanevski, J., Flores, R.O.R., Gabor, A. et al. Explainable multiview framework for dissecting spatial relationships from highly multiplexed data. Genome Biol 23, 97 (2022).
- [LIANA+ Documentation](https://liana-py.readthedocs.io/)
- [LIANA+ GitHub](https://github.com/saezlab/liana-py)

## Data Manager integration

Galaxy administrators install and run `data_manager_liana_resources`. It registers ligand-receptor, HCOP, and Metalinks files in the shared `liana_resources` Tool Data Table. Regular users do not run the Data Manager: the LIANA tool forms expose those server-side entries through `from_data_table`.

For `translate_resource`, the ligand-receptor source and the ortholog mapping source each support a cached mode. The ortholog selector also supports `builtin`, which invokes `li.resource.get_hcop_orthologs` for the selected target organism, and `history`, which accepts a user-supplied table.

### Tool Data Table and tests

The production table is declared in `tool_data_table_conf.xml.sample` and points to `tool-data/liana_resources.loc`; the `.loc.sample` file is only the installation template. For Planemo tests, `tool_data_table_conf.xml.test` points to `test-data/liana_resources.loc`, whose entries resolve local ligand-receptor, HCOP, and Metalinks fixtures through `${__HERE__}`. This mirrors the established IUC test-data-table pattern.

Repeated spatial-annulus, cell-type, HCOP-download, MuData, separator, and result-key parameters are defined in `macros.xml` so `cross_pcf`, `lric`, and plotting branches share identical defaults and validation.
