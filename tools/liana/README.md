# LIANA+ Galaxy Tool Suite (`liana 2.0.0`)

Galaxy tools for [LIANA+](https://liana-py.readthedocs.io/), a framework for single-cell, multi-sample, multi-modal, and spatial cell-cell communication (CCC) analysis.

---

## How the Tools Connect

```
┌────────────────────────────────────────────────────────────────────────────┐
│ 1. Liana Resource (resource.xml)                                           │
│    • Extract/translate L-R databases (.tsv)                                │
│    • Export MetalinksDB metabolite-receptor / enzyme / transporter tables   │
└───────────────────────────────┬────────────────────────────────────────────┘
                                │ (Custom L-R .tsv / MetalinksDB .tsv)
                                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ 2. Liana Utils (utils.xml)                                                 │
│    • Build spatial connectivity graphs: query_bandwidth → spatial_neighbors │
│    • Align/grid spatial slides: expand_coordinates / interpolate_adata      │
│    • Scale or merge modalities: zi_minmax / neg_to_zero / mdata_to_anndata  │
└───────────────────────────────┬────────────────────────────────────────────┘
                                │ (.h5ad / .h5mu with spatial connectivities)
                                ▼
┌─────────────────────────────────────────┬──────────────────────────────────┐
│ 3. Liana Methods (liana_methods.xml)    │ 4. Liana MISTy (liana_misty.xml) │
│    • Single-cell / multi-sample L-R     │    • Multi-view spatial learning │
│      (rank_aggregate, cellphonedb, ...) │      (genericMistyData,          │
│    • Spatial bivariate & inflow         │       lrMistyData, MistyData)    │
│    • Distance enrichment (lric,         │    • Outputs MISTy MuData (.h5mu)│
│      cross_pcf, get_lric_auc)           │      & metrics/interactions .tsv │
│    • DEA L-R (df_to_lr) &               │                                  │
│      metabolite (estimate_metalinks)    │                                  │
└───────────────┬─────────────────────────┴─────────────────┬────────────────┘
                │ (.h5ad with .uns[liana_res] or lrdata)    │
                ▼                                           │
┌─────────────────────────────────────────┐                 │
│ 5. Liana Multi-Sample (multi.xml)       │                 │
│    • Tensor-cell2cell (to_tensor_c2c)   │                 │
│    • MOFA views (lrs_to_views,          │                 │
│      adata_to_views, lrdata_to_mudata)  │                 │
│    • Spatial NMF (nmf) → extract scores │                 │
│      via Liana Utils (get_factor_scores,│                 │
│      get_variable_loadings)             │                 │
└───────────────┬─────────────────────────┘                 │
                │ (.h5ad / .h5mu)                           │ (.h5mu)
                ▼                                           ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ 6. Liana Plot (plot.xml)                                                   │
│    • L-R plots: dotplot, dotplot_by_sample, tileplot, circle               │
│    • Spatial & NMF plots: connectivity, feature_by_group, elbow            │
│    • Distance enrichment plots: annulus, lric_lineplot, lric_divergence    │
│    • MISTy plots: misty_target_metrics, misty_contributions,               │
│      misty_interactions                                                    │
└────────────────────────────────────────────────────────────────────────────┘
```

### Common End-to-End Workflows

1. **Single-Cell & Multi-Sample CCC**:
   - *(Optional)* **Liana Resource** (`get_hcop_orthologs` → `translate_resource`) to translate human L-R pairs to another species (`.tsv`).
   - **Liana Methods** (`rank_aggregate`, `cellphonedb`, `cellchat`, etc.) on an `AnnData` (`.h5ad`) or `MuData` (`.h5mu`) → outputs `.h5ad`/`.h5mu` with `.uns['liana_res']` (and `.tsv`).
   - *(Multi-sample)* **Liana Multi-Sample** (`to_tensor_c2c` for 4D tensor `.npy`/`.pkl`, or `lrs_to_views` / `adata_to_views` → `filter_view_markers` for `.h5mu` views).
   - **Liana Plot** (`dotplot`, `dotplot_by_sample`, `tileplot`, or `circle`) on the resulting `.h5ad`/`.h5mu`.

2. **Spatial Local Co-Expression, Inflow & NMF**:
   - **Liana Utils** (`query_bandwidth` → `spatial_neighbors`) on spatial `.h5ad` → outputs `.h5ad` with `.obsp['spatial_connectivities']` (inspectable with **Liana Plot** → `connectivity`).
   - **Liana Methods** (`bivariate` or `inflow` → `compute_global_specificity`) → outputs interaction `AnnData` (`lrdata`, `.h5ad`).
   - **Liana Plot** (`feature_by_group`) on `inflow` `.h5ad`, or **Liana Multi-Sample** (`lrdata_to_mudata` for `.h5mu` views, or `nmf` for spatial factor patterns followed by **Liana Plot** → `elbow` and **Liana Utils** → `get_factor_scores` / `get_variable_loadings`).

3. **Spatial Distance Enrichment (`lric` / `cross_pcf`)**:
   - **Liana Plot** (`annulus`) on spatial `.h5ad` to choose `radius_step` and `annulus_steps`.
   - **Liana Methods** (`lric` or `cross_pcf`) on `.h5ad` → stores `g(r)` curves in `.uns['lric']`, ranked via **Liana Methods** (`get_lric_auc`).
   - **Liana Plot** (`lric_lineplot` or `lric_divergence`) on the resulting `.h5ad`.

4. **Metabolite-Mediated CCC**:
   - **Liana Resource** (`get_metalinks` with `export_networks=True`) → outputs `metalinks_lr_resource.tsv`, `metalinks_pd_net.tsv`, and `metalinks_t_net.tsv`.
   - **Liana Methods** (`estimate_metalinks`) on `.h5ad` + the 3 MetalinksDB tables → outputs `.h5mu` (`metabolite` and `receptor` modalities), which can then be passed to **Liana Methods** (`rank_aggregate`, `bivariate`, etc.) or **Liana Utils** (`mdata_to_anndata`).

5. **Multi-View Spatial Modeling (MISTy)**:
   - *(Optional)* Prepare extra spatial views using **Liana Utils** (`obsm_to_adata`, `interpolate_adata`, `spatial_neighbors`).
   - **Liana MISTy** (`genericMistyData`, `lrMistyData`, or `MistyData`) on `.h5ad` or `.h5mu` → outputs MISTy `.h5mu` and `.tsv` tables.
   - **Liana Plot** (`misty_target_metrics`, `misty_contributions`, `misty_interactions`) on the MISTy `.h5mu`.

---

## Tool Reference (Functions, Inputs, and Outputs)

### 1. Liana Resource (`resource.xml` — `liana.rs`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `select_resource` | Extracts one of 17 built-in L-R resources (`consensus`, `mouseconsensus`, `cellphonedb`, etc.) | Built-in resource name | L-R table (`.tsv`) |
| `get_hcop_orthologs` | Downloads human-to-species gene orthologs from HGNC HCOP (19 species) | Target organism & min evidence | Ortholog mapping table (`.tsv`) |
| `translate_resource` | Translates ligand and receptor columns of an L-R resource across species | Built-in/custom L-R `.tsv` + ortholog `.tsv` | Translated L-R table (`.tsv`) |
| `translate_column` | Translates a single gene symbol column in any table | Input `.tsv` + ortholog `.tsv` | Translated table (`.tsv`) |
| `generate_lr_geneset` | Links L-R pairs (`ligand^receptor`) to pathway or disease gene sets | Built-in/custom L-R `.tsv` + network `.tsv` | L-R gene set table (`.tsv`) |
| `get_metalinks` | Queries MetalinksDB metabolite-protein interactions and exports networks for `estimate_metalinks` | Edge type & annotation filters | `metalinks.tsv` (+ optional `lr_resource.tsv`, `pd_net.tsv`, `t_net.tsv`) |

### 2. Liana Utils (`utils.xml` — `liana.pp` & `liana.ms`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `spatial_neighbors` | Computes spatial connectivity weights (or cross-modality rectangular weights) | `.h5ad` or `.h5mu` (+ optional reference `.h5ad`) | `.h5ad` or `.h5mu` with `.obsp['{key_added}_connectivities']` |
| `query_bandwidth` | Evaluates average neighbor counts across spatial bandwidths | `.h5ad` or `.h5mu` (+ optional reference `.h5ad`) | Bandwidth table (`.tsv`) + plot (`.png`/`.pdf`/`.svg`) |
| `spatial_pair_proximity` | Computes pairwise distances, proximity scores, and contact flags between cell types | `.h5ad` with `.obs[groupby]` & `.obsm[spatial_key]` | Pairwise proximity table (`.tsv`) |
| `expand_coordinates` | Lays out spatial coordinates of multiple samples onto a non-overlapping grid | `.h5ad` with `.obs[sample_key]` & `.obsm[spatial_key]` | `.h5ad` with grid `.obsm[spatial_key]` |
| `interpolate_adata` | Interpolates spatial features from a target `AnnData` onto reference coordinates | Target `.h5ad` + reference `.h5ad` | Interpolated `.h5ad` |
| `obsm_to_adata` | Extracts an `.obsm` matrix into a standalone `AnnData` | `.h5ad` or `.h5mu` | `.h5ad` |
| `mdata_to_anndata` | Concatenates two `MuData` modalities (`x_mod`, `y_mod`) into one `AnnData` | `.h5mu` | `.h5ad` |
| `zi_minmax` / `neg_to_zero` | Applies zero-inflated min-max scaling or clips values below cutoff to zero | `.h5ad` | `.h5ad` |
| `get_factor_scores` | Extracts observation factor scores from `.obsm` (NMF / MOFA) | `.h5ad` or `.h5mu` | Factor scores table (`.tsv`) |
| `get_variable_loadings` | Extracts feature loadings from `.varm` (NMF / MOFA) | `.h5ad` or `.h5mu` | Variable loadings table (`.tsv`) |

### 3. Liana Methods (`liana_methods.xml` — `liana.mt`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `rank_aggregate`, `cellphonedb`, `cellchat`, `natmi`, `connectome`, `logfc`, `singlecellsignalr`, `geometric_mean`, `scseqcomm` | Scores cell-type L-R interactions (single-sample, multi-sample `sample_key`, spatial proximity, or multi-modal `MuData`) | `.h5ad` or `.h5mu` (+ built-in/custom L-R resource) | `.h5ad` or `.h5mu` with `.uns[key_added]` (+ optional `.tsv`) |
| `bivariate` | Computes local (per-spot/cell) and global spatial bivariate L-R scores | `.h5ad` or `.h5mu` with `.obsp[connectivity_key]` | Interaction `lrdata` (`.h5ad`) (+ optional global `.tsv`) |
| `inflow` | Computes directed spatial inflow scores (`source^ligand^receptor`) per cell/spot | `.h5ad` with `.obsp[connectivity_key]` | Inflow `lrdata` (`.h5ad`) |
| `compute_global_specificity` | Computes permutation-based receiver cell-type specificity from Inflow scores | Inflow `.h5ad` | Updated `.h5ad` with `.uns['global_res']` (+ optional `.tsv`) |
| `lric` / `cross_pcf` | Computes distance-resolved pair correlation `g(r)` curves across concentric rings | `.h5ad` with `.obsm[spatial_key]` | `.h5ad` with `.uns[key_added]` (+ optional `.tsv`) |
| `get_lric_auc` | Summarizes `lric`/`cross_pcf` `g(r)` curves into AUC scores and peak radii | `.h5ad` with `.uns[uns_key]` | AUC summary table (`.tsv`) |
| `df_to_lr` | Combines per-cell-type DEA statistics with expression proportions into L-R scores | `.h5ad` + DEA table (`.tsv`) | Differential L-R table (`.tsv`) |
| `estimate_metalinks` | Estimates metabolite abundance from enzyme/transporter expression and pairs with receptors | `.h5ad` + MetalinksDB `resource`, `pd_net` (& optional `t_net`) `.tsv` | `.h5mu` (`metabolite` & `receptor` modalities) |

### 4. Liana MISTy (`liana_misty.xml` — `liana.mt`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `genericMistyData` | Builds `intra`, `juxta`, and `para` views (plus optional separate `extra` modality) and fits MISTy | `.h5ad` (+ optional `extra` `.h5ad`) | MISTy `.h5mu` (+ optional `target_metrics.tsv` & `interactions.tsv`) |
| `lrMistyData` | Builds receptor `intra` and spatially weighted ligand `extra` views and fits MISTy | `.h5ad` (+ built-in/custom L-R resource) | MISTy `.h5mu` (+ optional `target_metrics.tsv` & `interactions.tsv`) |
| `MistyData` | Fits MISTy on custom extra-view `AnnData` files or an existing `MuData` | `.h5ad` + custom view `.h5ad` files, or `.h5mu` | MISTy `.h5mu` (+ optional `target_metrics.tsv` & `interactions.tsv`) |

### 5. Liana Multi-Sample (`multi.xml` — `liana.ms`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `to_tensor_c2c` | Reshapes multi-sample L-R scores into a 4D tensor (`samples × senders × receivers × LRs`) | `.h5ad` with by-sample `.uns[uns_key]` | 4D tensor (`.npy` or `.pkl`) |
| `lrs_to_views` | Converts multi-sample L-R scores into `source&target` views for MOFAtalk | `.h5ad` with by-sample `.uns[uns_key]` | `.h5mu` |
| `adata_to_views` | Aggregates single-cell data into pseudobulk views per cell type for MOFAcellular | `.h5ad` with `.obs[groupby]` & `.obs[sample_key]` | `.h5mu` |
| `filter_view_markers` | Filters background marker genes of other cell types from `adata_to_views` modalities | `.h5mu` + marker gene table (`.tsv`) | Filtered `.h5mu` |
| `lrdata_to_mudata` | Splits a spatial Inflow `AnnData` into one `MuData` view per sender cell type | Inflow `.h5ad` | `.h5mu` |
| `nmf` | Runs NMF on an `AnnData` matrix (fixed `n_components` or automatic elbow rank estimation) | `.h5ad` | `.h5ad` with `.obsm['NMF_W']`, `.varm['NMF_H']`, `.uns['nmf_errors']` |

### 6. Liana Plot (`plot.xml` — `liana.pl`)

| Function | What it does | Input | Output |
| :--- | :--- | :--- | :--- |
| `dotplot` / `dotplot_by_sample` | Dotplot of L-R interactions (single-sample or faceted by sample) | `.h5ad` or `.h5mu` with `.uns[uns_key]` | `.png`, `.pdf`, or `.svg` |
| `tileplot` | Side-by-side ligand and receptor tiles colored/labeled by column suffix | `.h5ad` or `.h5mu` with `.uns[uns_key]` | `.png`, `.pdf`, or `.svg` |
| `circle` | Circular network plot of interaction counts or mean scores between cell types | `.h5ad` or `.h5mu` with `.uns[uns_key]` | `.png`, `.pdf`, or `.svg` |
| `connectivity` | Spatial plot of neighbor connectivity weights for a selected cell/spot | `.h5ad` or `.h5mu` with `.obsp[connectivity_key]` | `.png`, `.pdf`, or `.svg` |
| `elbow` | Plot of NMF reconstruction error per rank (`k`) and estimated elbow rank | `.h5ad` with `.uns['nmf_errors']` | `.png`, `.pdf`, or `.svg` |
| `feature_by_group` | Spatial scatter of an Inflow feature (`source^ligand^receptor`) across receiver groups | Inflow `.h5ad` | `.png`, `.pdf`, or `.svg` |
| `annulus` | Overlays concentric distance rings on spatial coordinates to select `lric`/`cross_pcf` radii | `.h5ad` with `.obsm[spatial_key]` | `.png`, `.pdf`, or `.svg` |
| `lric_lineplot` / `lric_divergence` | Plots `g(r)` distance enrichment curves or compares two `log2 g(r)` curves | `.h5ad` or `.h5mu` with `.uns[uns_key]` | `.png`, `.pdf`, or `.svg` |
| `misty_target_metrics` / `misty_contributions` / `misty_interactions` | Plots MISTy R² metrics, view contributions, or predictor-target importances | MISTy `.h5mu` | `.png`, `.pdf`, or `.svg` |

---

## References

- Dimitrov, D., Schäfer, P.S.L., Farr, E. et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nat Cell Biol* 26, 1613–1622 (2024). https://doi.org/10.1038/s41556-024-01469-w
- Dimitrov, D., Türei, D., Garrido-Rodriguez, M. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun* 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0
- Tanevski, J., Flores, R.O.R., Gabor, A. et al. Explainable multiview framework for dissecting spatial relationships from highly multiplexed data. *Genome Biol* 23, 97 (2022). https://doi.org/10.1186/s13059-022-02663-5
