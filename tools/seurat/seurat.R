options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dplyr)
    library(optparse)
    library(rmarkdown)
    library(ggplot2)
})

min_cells <- 3
min_genes <- 200
low_thresholds <- 1
high_thresholds <- 20000000
x_low_cutoff <- 0.0125
x_high_cutoff <- 3
y_cutoff <- 0.5
numPCs <- 10
cells_use <- 500
resolution <- 0.6
min_pct <- 0.25
logfc_threshold <- 0.25


# *Read in Data, generate inital Seurat object*

counts <- read.delim("~/Downloads/deng_small.tab.gz", row.names=1)
seuset <- CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_genes)


# Raw data vizualization
VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"), axis="v")
FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# Filter and normalize for UMI counts

seuset <- subset(seuset, subset = `nCount_RNA` > low_thresholds & `nCount_RNA` < high_thresholds)
seuset <- NormalizeData(seuset, normalizeation.method = "LogNormalize", scale.factor = 10000)

# Variable Genes

seuset <- FindVariableFeatures(seuset,selection.method = "mvp",
                               x.low.cutoff = x_low_cutoff,
                               x.high.cutoff = x_high_cutoff,
                               y.cutoff = y_cutoff
                               )

VariableFeaturePlot(seuset, cols = c("black", "red"), selection.method = "disp")


seuset <- ScaleData(object = seuset, vars.to.regress = "nCount_RNA")

# PCA vizualization
seuset <- RunPCA(seuset, npcs=10)
VizDimLoadings(seuset, dims = 1:2)
PCAPlot(seuset, dims = c(1,2))
DimHeatmap(seuset, dims=1:10, nfeatures=30, reduction="pca")
seuset <- JackStraw(seuset, dims=10, reduction = "pca", num.replicate = 100)
seuset <- ScoreJackStraw(seuset, dims = 1:10)
JackStrawPlot(seuset, dims = 1:10)
ElbowPlot(seuset, ndims = 20, reduction = "pca")

# tSNE
seuset <- FindNeighbors(object = seuset)
seuset <- FindClusters(object = seuset)
seuset <- RunTSNE(seuset, dims = 1:10, resolution =0.6)
TSNEPlot(seuset)

# Marker Genes
markers <- FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(seuset, features = top10$gene)

rmarkdown::render('/Users/alexanderostrovsky/Desktop/galaxy/tools/myTools/seurat/Seurat.R')
