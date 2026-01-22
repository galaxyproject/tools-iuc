LIANA+
======

Galaxy tools for [LIANA+](https://liana-py.readthedocs.io/), a framework for cell-cell communication analysis.

## Tools

### 1. Single Cell Methods (`single_cell_copy.xml`)

Ligand-receptor inference methods for single-cell data.

Method | Description
--- | ---
`cellchat` | CellChat ligand-receptor method
`cellphonedb` | CellPhoneDB ligand-receptor method
`connectome` | Connectome ligand-receptor method
`logfc` | Log fold change ligand-receptor method
`natmi` | NATMI ligand-receptor method
`singlecellsignalr` | SingleCellSignalR ligand-receptor method
`geometric_mean` | Geometric mean ligand-receptor method
`rank_aggregate` | Aggregate rankings from multiple methods

### 2. Spatial Methods (`spatial_copy.xml`)

Spatial ligand-receptor analysis and MISTy modelling.

Method | Description
--- | ---
`bivariate` | Bivariate local spatial metrics for ligand-receptor analysis
`genericMistyData` | Generic MISTy multi-view modelling
`lrMistyData` | Ligand-receptor MISTy modelling

#### Bivariate Metrics

**Local metrics:**

Metric | Description
--- | ---
`cosine` | Weighted cosine similarity
`jaccard` | Weighted Jaccard similarity
`pearson` | Weighted Pearson correlation
`spearman` | Weighted Spearman correlation
`masked_spearman` | Masked and weighted Spearman correlation
`product` | Simple weighted product
`norm_product` | Normalized weighted product
`morans` | Moran's R

**Global metrics:**

Metric | Description
--- | ---
`morans` | Moran's I (global spatial autocorrelation)
`lee` | Lee's L (bivariate spatial association)

#### MISTy Models

Model | Description
--- | ---
`RandomForestModel` | Random Forest (uses out-of-bag predictions)
`LinearModel` | Linear regression (uses cross-validation)
`RobustLinearModel` | Robust linear regression

#### Kernel Functions

Kernel | Description
--- | ---
`misty_rbf` | MISTy RBF kernel (Gaussian derivative)
`gaussian` | Gaussian kernel
`exponential` | Exponential kernel
`linear` | Linear kernel

### 3. Available Resources

Built-in ligand-receptor resources from `li.resource.show_resources()`:

- `consensus` (default)
- `baccin2019`
- `cellcall`
- `cellchatdb`
- `cellinker`
- `cellphonedb`
- `celltalkdb`
- `connectomedb2020`
- `embrace`
- `guide2pharma`
- `hpmr`
- `icellnet`
- `italk`
- `kirouac2010`
- `lrdb`
- `mouseconsensus`
- `ramilowski2015`

## References

- Dimitrov, D., Türei, D., Garrido-Rodriguez, M. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0
- [LIANA+ Documentation](https://liana-py.readthedocs.io/)
- [LIANA+ GitHub](https://github.com/saezlab/liana-py)
