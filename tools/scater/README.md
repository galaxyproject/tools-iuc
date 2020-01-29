# Wrappers for Scater

This code wraps a number of [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) functions as Galaxy wrappers. Briefly, the `scater-create-qcmetric-ready-sce` tool takes a sample gene expression matrix (usually read-counts) and a cell annotation file, creates a [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object and runs scater's `calculateQCMetrics` function (using other supplied files such as ERCC's and mitochondrial gene features).
Various filter scripts are provided, along with some plotting functions for QC.


## Typical workflow

1. Read in data with `scater-create-qcmetric-ready-sce`.
2. Visualise it.\
   Take a look at the distribution of library sizes, expressed features and mitochondrial genes with `scater-plot-dist-scatter`.
   Then look at the distibution of genes across cells with `scater-plot-exprs-freq`.
3. Guided by the plots, filter the data with `scater-filter`.\
   You can either manually filter with user-defined parameters or use PCA to automatically removes outliers.
4. Visualise data again to see how the filtering performed using `scater-plot-dist-scatter`.\
   Decide if you're happy with the data. If not, try increasing or decreasing the filtering parameters.
5. Normalise data with `scater-normalize`.
6. Investigate other confounding factors.\
   Plot the data (using PCA) and display various annotated properties of the cells using `scater-plot-pca`.

## Command-line usage

The scripts require the installation of scater and few other R/BioConductor packages. An easy way to install them is to create a [conda](https://conda.io/) environment using the `environment.yml` file distributed together with these wrappers:

```
conda env create -f environment.yml
conda activate scater
```

For help with any of the following scripts, run:
 `<script-name> --help`

---

`scater-create-qcmetric-ready-sce.R`
Takes an expression matrix (usually read-counts) of samples (columns) and gene/transcript features (rows), along with other annotation information, such as cell metadata, control genes (mitochondrail genes, ERCC's), creates a [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object and runs scater's `calculateQCMetrics`. Save the resulting SingleCellExperiment object in Loom format.


```
./scater-create-qcmetric-ready-sce.R -a test-data/counts.txt -c test-data/annotation.txt -f test-data/mt_controls.txt  -o test-data/scater_qcready.loom
```

---

`scater-plot-dist-scatter.R`
Takes SingleCellExperiment object (from Loom file) and plots a panel of read and feature graphs, including the distribution of library sizes, distribution of feature counts, a scatterplot of reads vs features, and % of mitochondrial genes in library.

```
./scater-plot-dist-scatter.R -i test-data/scater_qcready.loom -o test-data/scater_reads_genes_dist.pdf
```

---

`scater-plot-exprs-freq.R`
Plots mean expression vs % of expressing cells and provides information as to the number of genes expressed in 50% and 25% of cells.

---

`scater-pca-filter.R`
Takes SingleCellExperiment object (from Loom file) and automatically removes outliers from data using PCA. Save the filtered SingleCellExperiment object in Loom format.

```
./scater-pca-filter.R -i test-data/scater_qcready.loom -o test-data/scater_pca_filtered.loom
```

---

`scater-manual-filter.R`
Takes SingleCellExperiment object (from Loom file) and filters data using user-provided parameters. Save the filtered SingleCellExperiment object in Loom format.

```
./scater-manual-filter.R -i test-data/scater_qcready.loom -l 10000 -d 4 -m 33 -o test-data/scater_manual_filtered.loom
```

---

`scater-normalize.R`
Compute log-normalized expression values from count data in a SingleCellExperiment object, using the size factors stored in the object. Save the normalised SingleCellExperiment object in Loom format.

```
./scater-normalize.R -i test-data/scater_manual_filtered.loom -o test-data/scater_man_filtered_normalised.loom
```

---

`scater-plot-pca.R`
PCA plot of a normalised SingleCellExperiment object (produced with `scater-normalize.R`). The options `-c`, `-p`, and `-s` all refer to cell annotation features. These are the column headers of the `-c` option used in `scater-create-qcmetric-ready-sce.R`.

```
./scater-plot-pca.R -i test-data/scater_man_filtered_normalised.loom -c Treatment -p Mutation_Status -o test-data/scater_pca_plot.pdf
```
