#' ---
#' title: "Seurat Analysis"
#' author: "Performed using Galaxy"
#' params:
#'     counts: ""
#'     min_cells: ""
#'     min_genes: ""
#'     low_thresholds: ""
#'     high_thresholds: ""
#'     numPCs: ""
#'     cells_use: ""
#'     resolution: ""
#'     perplexity: ""
#'     min_pct: ""
#'     logfc_threshold: ""
#'     showcode: ""
#'     warn: ""
#'     varstate: ""
#'     vlnfeat: ""
#'     featplot: ""
#'     PCplots: ""
#'     tsne: ""
#'     heatmaps: ""
#' ---

# nolint start
#+ echo=F, warning = F, message=F
options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr()); q("no", 1, F)
})
showcode <- as.logical(params$showcode)
warn <-  as.logical(params$warn)
varstate <- as.logical(params$varstate)
vlnfeat <- as.logical(params$vlnfeat)
featplot <- as.logical(params$featplot)
pc_plots <- as.logical(params$PCplots)
tsne <- as.logical(params$tsne)
heatmaps <- as.logical(params$heatmaps)

# we need that to not crash Galaxy with an UTF-8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

#+ echo = F, warning = `warn`, include =`varstate`
min_cells <- as.integer(params$min_cells)
min_genes <- as.integer(params$min_genes)
low_thresholds <- as.integer(params$low_thresholds)
high_thresholds <- as.integer(params$high_thresholds)
num_pcs <- as.integer(params$numPCs)
cells_use <- as.integer(params$cells_use)
resolution <- as.double(params$resolution)
perplexity <- as.integer(params$perplexity)
min_pct <- as.double(params$min_pct)
logfc_threshold <- as.double(params$logfc_thresh)
print(paste0("Minimum cells: ", min_cells))
print(paste0("Minimum features: ", min_genes))
print(paste0("Umi low threshold: ", low_thresholds))
print(paste0("Umi high threshold: ", high_thresholds))
print(paste0("Number of principal components: ", num_pcs))
print(paste0("Resolution: ", resolution))
print(paste0("Perplexity: ", perplexity))
print(paste0("Minimum percent of cells", min_pct))
print(paste0("Logfold change threshold", logfc_threshold))

#+ echo = FALSE
if (showcode == TRUE) print("Read in data, generate inital Seurat object")
#+ echo = `showcode`, warning = `warn`, message = F
counts <- read.delim(params$counts, row.names = 1)
seuset <- Seurat::CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_genes)

#+ echo = FALSE
if (showcode == TRUE && vlnfeat == TRUE) print("Raw data vizualization")
#+ echo = `showcode`, warning = `warn`, include=`vlnfeat`
Seurat::VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"))
Seurat::FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#+ echo = FALSE
if (showcode == TRUE) print("Filter and normalize for UMI counts")
#+ echo = `showcode`, warning = `warn`
seuset <- subset(seuset, subset = `nCount_RNA` > low_thresholds & `nCount_RNA` < high_thresholds)
seuset <- Seurat::NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)

#+ echo = FALSE
if (showcode == TRUE && featplot == TRUE) print("Variable Genes")
#+ echo = `showcode`, warning = `warn`, include = `featplot`
seuset <- Seurat::FindVariableFeatures(object = seuset, selection.method = "mvp")
Seurat::VariableFeaturePlot(seuset, cols = c("black", "red"), selection.method = "disp")
seuset <- Seurat::ScaleData(object = seuset, vars.to.regress = "nCount_RNA")

#+ echo = FALSE
if (showcode == TRUE && pc_plots == TRUE) print("PCA Visualization")
#+ echo = `showcode`, warning = `warn`, include = `pc_plots`
seuset <- Seurat::RunPCA(seuset, npcs = num_pcs)
Seurat::VizDimLoadings(seuset, dims = 1:2)
Seurat::DimPlot(seuset, dims = c(1, 2), reduction = "pca")
Seurat::DimHeatmap(seuset, dims = 1:num_pcs, nfeatures = 30, reduction = "pca")
seuset <- Seurat::JackStraw(seuset, dims = num_pcs, reduction = "pca", num.replicate = 100)
seuset <- Seurat::ScoreJackStraw(seuset, dims = 1:num_pcs)
Seurat::JackStrawPlot(seuset, dims = 1:num_pcs)
Seurat::ElbowPlot(seuset, ndims = num_pcs, reduction = "pca")

#+ echo = FALSE
if (showcode == TRUE && tsne == TRUE) print("tSNE")
#+ echo = `showcode`, warning = `warn`, include = `tsne`
seuset <- Seurat::FindNeighbors(object = seuset)
seuset <- Seurat::FindClusters(object = seuset)
if (perplexity == -1) {
    seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution);
} else {
    seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution, perplexity = perplexity);
}
Seurat::DimPlot(seuset, reduction = "tsne")

#+ echo = FALSE
if (showcode == TRUE && heatmaps == TRUE) print("Marker Genes")
#+ echo = `showcode`, warning = `warn`, include = `heatmaps`
markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)
top10 <- dplyr::group_by(markers, cluster)
top10 <- dplyr::top_n(top10, 10, avg_log2FC)
Seurat::DoHeatmap(seuset, features = top10$gene)
# nolint end
