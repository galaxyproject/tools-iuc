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
#'     min_pct: ""
#'     logfc_threshold: ""
#'     showcode: ""
#'     warn: ""
#'     varstate: ""
#' ---

#+ echo=F, warning = F, message=F
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
showcode <- as.logical(params$showcode)
warn <-  as.logical(params$warn)
varstate <- as.logical(params$varstate)
# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")


#+ echo = F, warning = `warn`, include =`varstate`
min_cells <- as.integer(params$min_cells)
min_genes <- as.integer(params$min_genes)
low_thresholds <- as.integer(params$low_thresholds)
high_thresholds <- as.integer(params$high_thresholds)
numPCs <- as.integer(params$num_PCs)
cells_use <- as.integer(params$cells_use)
resolution <- as.double(params$resolution)
min_pct <- as.double(params$min_pct)
logfc_threshold <- as.double(params$logfc_thresh)
print(paste0("Minimum cells: ", min_cells))
print(paste0("Minimum features: ", min_genes))
print(paste0("Umi low threshold: ", low_thresholds))
print(paste0("Umi high threshold: ", high_thresholds))


#+ echo = FALSE
if(showcode == TRUE){print("Read in data, generate inital Seurat object")}
#+ echo = `showcode`, warning = `warn`, message = F
counts <- read.delim(params$counts, row.names=1)
seuset <- Seurat::CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_genes)

#+ echo = FALSE
if(showcode == TRUE){print("Raw data vizualization")}
#+ echo = `showcode`, warning = `warn`
Seurat::VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"), axis="v")
Seurat::FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#+ echo = FALSE
if(showcode == TRUE){print("Filter and normalize for UMI counts")}
#+ echo = `showcode`, warning = `warn`
seuset <- subset(seuset, subset = `nCount_RNA` > low_thresholds & `nCount_RNA` < high_thresholds)
seuset <- Seurat::NormalizeData(seuset, normalizeation.method = "LogNormalize", scale.factor = 10000)

#+ echo = FALSE
if(showcode == TRUE){print("Variable Genes")}
#+ echo = `showcode`, warning = `warn`
seuset <- Seurat::FindVariableFeatures(object = seuset, selection.method = "mvp")
Seurat::VariableFeaturePlot(seuset, cols = c("black", "red"), selection.method = "disp")
seuset <- Seurat::ScaleData(object = seuset, vars.to.regress = "nCount_RNA")

#+ echo = FALSE
if(showcode == TRUE){print("PCA Visualization")}
#+ echo = `showcode`, warning = `warn`
seuset <- Seurat::RunPCA(seuset, npcs=10)
Seurat::VizDimLoadings(seuset, dims = 1:2)
Seurat::DimPlot(seuset, dims = c(1,2), reduction="pca")
Seurat::DimHeatmap(seuset, dims=1:10, nfeatures=30, reduction="pca")
seuset <- Seurat::JackStraw(seuset, dims=10, reduction = "pca", num.replicate = 100)
seuset <- Seurat::ScoreJackStraw(seuset, dims = 1:10)
Seurat::JackStrawPlot(seuset, dims = 1:10)
Seurat::ElbowPlot(seuset, ndims = 20, reduction = "pca")

#+ echo = FALSE
if(showcode == TRUE){print("tSNE")}
#+ echo = `showcode`, warning = `warn`
seuset <- Seurat::FindNeighbors(object = seuset)
seuset <- Seurat::FindClusters(object = seuset)
seuset <- Seurat::RunTSNE(seuset, dims = 1:10, resolution =0.6)
Seurat::DimPlot(seuset, reduction="tsne")

#+ echo = FALSE
if(showcode == TRUE){print("Marker Genes")}
#+ echo = `showcode`, warning = `warn`
markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)
top10 <- dplyr::group_by(markers, cluster)
top10 <- dplyr::top_n(top10, 10, avg_logFC)
Seurat::DoHeatmap(seuset, features = top10$gene)
