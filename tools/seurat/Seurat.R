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
#'     resolution: ""
#'     perplexity: ""
#'     min_pct: ""
#'     logfc_threshold: ""
#'     start_step: ""
#'     end_step: ""
#'     showcode: ""
#'     warn: ""
#'     varstate: ""
#'     vlnfeat: ""
#'     featplot: ""
#'     PCplots: ""
#'     tsne: ""
#'     heatmaps: ""
#'     norm_out: ""
#'     variable_out: ""
#'     pca_out : ""
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
start_step <- as.integer(params$start_step)
end_step <- as.integer(params$end_step)
norm_out <- as.logical(params$norm_out)
# we need that to not crash Galaxy with an UTF-8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

#+ echo = F, warning = `warn`, include =`varstate`
min_cells <- as.integer(params$min_cells)
min_genes <- as.integer(params$min_genes)
low_thresholds <- as.integer(params$low_thresholds)
high_thresholds <- as.integer(params$high_thresholds)

if (end_step >= 2) {
    variable_out <- as.logical(params$variable_out)
}


if (end_step > 3) {
    num_pcs <- as.integer(params$numPCs)
    pca_out <- as.logical(params$pca_out)
}
if (end_step > 4) {
    perplexity <- as.integer(params$perplexity)
    resolution <- as.double(params$resolution)
    clusters_out <- as.logical(params$clusters_out)
}
if (end_step > 5) {
    min_pct <- as.double(params$min_pct)
    logfc_threshold <- as.double(params$logfc_thresh)
    markers_out <- as.logical(params$markers_out)
}

# min_cells <- as.integer(params$min_cells)
# min_genes <- as.integer(params$min_genes)
# low_thresholds <- as.integer(params$low_thresholds)
# high_thresholds <- as.integer(params$high_thresholds)
# num_pcs <- as.integer(params$numPCs)
# resolution <- as.double(params$resolution)
# perplexity <- as.integer(params$perplexity)
# min_pct <- as.double(params$min_pct)
# logfc_threshold <- as.double(params$logfc_thresh)
# print(paste0("Minimum cells: ", min_cells))
# print(paste0("Minimum features: ", min_genes))
# print(paste0("Umi low threshold: ", low_thresholds))
# print(paste0("Umi high threshold: ", high_thresholds))
# print(paste0("Number of principal components: ", num_pcs))
# print(paste0("Resolution: ", resolution))
# print(paste0("Perplexity: ", perplexity))
# print(paste0("Minimum percent of cells", min_pct))
# print(paste0("Logfold change threshold", logfc_threshold))


if (start_step == 0) {
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

}

if (start_step <= 1) {
    #+ echo = FALSE
    if (showcode == TRUE) print("Filter and normalize for UMI counts")
    #+ echo = `showcode`, warning = `warn`
    seuset <- subset(seuset, subset = `nCount_RNA` > low_thresholds & `nCount_RNA` < high_thresholds)
    seuset <- Seurat::NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
    if (norm_out == TRUE) {
         saveRDS(seuset, "norm_out.rds")

    }
}

if (end_step >= 2) {
    #+ echo = FALSE
    if (showcode == TRUE && featplot == TRUE) print("Variable Genes")
    #+ echo = `showcode`, warning = `warn`, include = `featplot`
    seuset <- Seurat::FindVariableFeatures(object = seuset, selection.method = "mvp")
    Seurat::VariableFeaturePlot(seuset, cols = c("black", "red"), selection.method = "disp")
    seuset <- Seurat::ScaleData(object = seuset, vars.to.regress = "nCount_RNA")
    if (variable_out == TRUE) {
        saveRDS(seuset, "var_out.rds")
    }
}

if (end_step >= 3) {
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
    if (pca_out == TRUE) {
        # pcadata <- Embeddings(seuset, reduction="pca")
        # fwrite(x = scaledata, sep="\t", file = "pcadata.tsv")
        Seurat::saveRDS(seuset, "pca_out.rds")
    }
}

if (end_step >= 4) {
#     #+ echo = FALSE
#     if (showcode == TRUE && tsne == TRUE) print("tSNE")
#     #+ echo = `showcode`, warning = `warn`, include = `tsne`
#     seuset <- Seurat::FindNeighbors(object = seuset)
#     seuset <- Seurat::FindClusters(object = seuset)
#     if (perplexity == -1) {
#         seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution);
#     } else {
#         seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution, perplexity = perplexity);
#     }
#     Seurat::DimPlot(seuset, reduction = "tsne")
#     if (tsne_out == TRUE) {
#         tsnedata <- Embeddings(seuset, reduction="tsne")
#         fwrite(x = tsnedata, sep="\t", file = "tsnedata.tsv")
#     }
# }
    #+ echo = FALSE
    if (showcode == TRUE && tsne == TRUE) print("tSNE and UMAP")
    #+ echo = `showcode`, warning = `warn`, include = `tsne`
    seuset <- Seurat::FindNeighbors(object = seuset)
    seuset <- Seurat::FindClusters(object = seuset)
    if (perplexity == -1) {
        seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution);
    } else {
        seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution, perplexity = perplexity);
    }
    Seurat::DimPlot(seuset, reduction = "tsne")
    seuset <- Seurat::RunUMAP(seuset, dims = 1:num_pcs)
    Seurat::DimPlot(seuset, reduction = "umap")
    if (clustersout == TRUE) {
        tsnedata <- Seurat::Embeddings(seuset, reduction="tsne")
        Seurat::saveRDS(seuset, "tsne_out.rds")
        umapdata <- Seurat::Embeddings(seuset, reduction="umap")
        Seurat::saveRDS(seuset, "umap_out.rds")
    }
}


if (end_step == 5) {
    #+ echo = FALSE
    if (showcode == TRUE && heatmaps == TRUE) print("Marker Genes")
    #+ echo = `showcode`, warning = `warn`, include = `heatmaps`
    markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)
    top10 <- dplyr::group_by(markers, cluster)
    top10 <- dplyr::top_n(top10, 10, avg_log2FC)
    Seurat::DoHeatmap(seuset, features = top10$gene)
    if (marker_out == TRUE) {
        Seurat::saveRDS(seuset, "markers_out.rds")
    }
}
# nolint end
