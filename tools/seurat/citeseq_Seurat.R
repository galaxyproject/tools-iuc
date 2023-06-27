#' ---
#' title: "Seurat Cite-seq Analysis"
#' author: "Performed using Galaxy"
#' params:
#'     rna: ""
#'     prot: ""
#'     min_cells: ""
#'     min_genes: ""
#'     low_thresholds: ""
#'     high_thresholds: ""
#'     numPCs: ""
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
#'     nmds: ""
#'     heatmaps: ""
#'     norm_out: ""
#'     variable_out: ""
#'     pca_out : ""
#'     clusters_out: ""
#'     markers_out: ""
#'     cite_markers: ""
#'     comparison: ""
#'     feat_comp: ""
#'     marker_compare: ""
#'     top_x: ""
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
nmds <- as.logical(params$nmds)
heatmaps <- as.logical(params$heatmaps)
end_step <- as.integer(params$end_step)
norm_out <- as.logical(params$norm_out)
comparison <- as.logical(params$comparison)
feature <- trimws(unlist(strsplit(as.character(params$feat_comp), ",")))
marker_compare <- as.logical(params$marker_compare)
top_x <- as.integer(params$top_x)


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

# we need that to not crash Galaxy with an UTF-8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

#+ echo = F, warning = `warn`, include =`varstate`
min_cells <- as.integer(params$min_cells)
min_genes <- as.integer(params$min_genes)
print(paste0("Minimum cells: ", min_cells))
print(paste0("Minimum features: ", min_genes))
low_thresholds <- as.integer(params$low_thresholds)
high_thresholds <- as.integer(params$high_thresholds)
print(paste0("Umi low threshold: ", low_thresholds))
print(paste0("Umi high threshold: ", high_thresholds))
variable_out <- as.logical(params$variable_out)
num_pcs <- as.integer(params$numPCs)
print(paste0("Number of principal components: ", num_pcs))
pca_out <- as.logical(params$pca_out)


if (params$perplexity == "") {
    perplexity <- -1
    print(paste0("Perplexity: ", perplexity))
} else { 
    perplexity <- as.integer(params$perplexity)
    print(paste0("Perplexity: ", perplexity))
}

resolution <- as.double(params$resolution)
print(paste0("Resolution: ", resolution))
clusters_out <- as.logical(params$clusters_out)

min_pct <- as.double(params$min_pct)
logfc_threshold <- as.double(params$logfc_thresh)
print(paste0("Minimum percent of cells", min_pct))
print(paste0("Logfold change threshold", logfc_threshold))
markers_out <- as.logical(params$markers_out)


#+ echo = FALSE
if (showcode == TRUE) print("Read in data, generate inital Seurat object")
#+ echo = `showcode`, warning = `warn`, message = F
rna <- read.delim(params$rna, row.names = 1)
rna <- Seurat::CollapseSpeciesExpressionMatrix(rna)
protein <- read.delim(params$prot, row.names = 1)
tryCatch(all.equal(colnames(rna), colnames(protein)), error = "Columns do not match in input files")
seuset <- Seurat::CreateSeuratObject(counts = rna, min.cells = min_cells, min.features = min_genes)

if (showcode == TRUE) print("asdf")
#+ echo = `showcode`, warning = `warn`, message = F
prot_obj <- Seurat::CreateAssayObject(counts = protein)

if (showcode == TRUE) print("qwer")
#+ echo = `showcode`, warning = `warn`, message = F
seuset[["ADT"]] <- prot_obj

if (showcode == TRUE) print("zxcv")
#+ echo = `showcode`, warning = `warn`, message = F
Seurat::DefaultAssay(seuset) <- "RNA"

if (showcode == TRUE && vlnfeat == TRUE) print("Raw data vizualization")
#+ echo = `showcode`, warning = `warn`, include=`vlnfeat`
if (vlnfeat == TRUE){
    print(Seurat::VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA")))
    print(Seurat::FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
}

if (showcode == TRUE) print("Filter and normalize for UMI counts")
#+ echo = `showcode`, warning = `warn`
seuset <- subset(seuset, subset = `nCount_RNA` > low_thresholds & `nCount_RNA` < high_thresholds)
seuset <- Seurat::NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
if (norm_out == TRUE) {
        saveRDS(seuset, "norm_out.rds")
}


if (showcode == TRUE && featplot == TRUE) print("Variable Genes")
#+ echo = `showcode`, warning = `warn`, include = `featplot`
seuset <- Seurat::FindVariableFeatures(object = seuset, selection.method = "mvp")
if (featplot == TRUE) {
    print(Seurat::VariableFeaturePlot(seuset, cols = c("black", "red"), selection.method = "disp"))
}
seuset <- Seurat::ScaleData(object = seuset, vars.to.regress = "nCount_RNA")
if (variable_out == TRUE) {
    saveRDS(seuset, "var_out.rds")
}



if (showcode == TRUE && pc_plots == TRUE) print("PCA Visualization")
#+ echo = `showcode`, warning = `warn`, include = `pc_plots`
seuset <- Seurat::RunPCA(seuset, npcs = num_pcs)
seuset <- Seurat::JackStraw(seuset, dims = num_pcs, reduction = "pca", num.replicate = 100)
seuset <- Seurat::ScoreJackStraw(seuset, dims = 1:num_pcs)
if (pc_plots == TRUE) {
    print(Seurat::VizDimLoadings(seuset, dims = 1:2))
    print(Seurat::DimPlot(seuset, dims = c(1, 2), reduction = "pca"))
    print(Seurat::DimHeatmap(seuset, dims = 1:num_pcs, nfeatures = 10, reduction = "pca"))
    print(Seurat::JackStrawPlot(seuset, dims = 1:num_pcs))
    print(Seurat::ElbowPlot(seuset, ndims = num_pcs, reduction = "pca"))
}
if (pca_out == TRUE) {
    saveRDS(seuset, "pca_out.rds")
}



if (showcode == TRUE && nmds == TRUE) print("tSNE and UMAP")
#+ echo = `showcode`, warning = `warn`, include = `nmds`
seuset <- Seurat::FindNeighbors(object = seuset)
seuset <- Seurat::FindClusters(object = seuset)
if (perplexity == -1) {
    seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution, check_duplicates = FALSE)
} else {
    seuset <- Seurat::RunTSNE(seuset, dims = 1:num_pcs, resolution = resolution, perplexity = perplexity, check_duplicates = FALSE)
}
if (nmds == TRUE) {
    print(Seurat::DimPlot(seuset, reduction = "tsne"))
}
seuset <- Seurat::RunUMAP(seuset, dims = 1:num_pcs)
if (nmds == TRUE) {
        print(Seurat::DimPlot(seuset, reduction = "umap"))
}
if (clusters_out == TRUE) {
    tsnedata <- Seurat::Embeddings(seuset, reduction="tsne")
    saveRDS(seuset, "tsne_out.rds")
    umapdata <- Seurat::Embeddings(seuset, reduction="umap")
    saveRDS(seuset, "umap_out.rds")
}

if (showcode == TRUE && heatmaps == TRUE) print("Marker Genes")
#+ echo = `showcode`, warning = `warn`, include = `heatmaps`
markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)
top10 <- dplyr::group_by(markers, cluster)
top10 <- dplyr::top_n(top10, n = 10, wt = avg_log2FC)
print(top10)
if (heatmaps == TRUE) {
    print(Seurat::DoHeatmap(seuset, features = top10$gene))
}
if (markers_out == TRUE) {
    saveRDS(seuset, "markers_out.rds")
    data.table::fwrite(x = markers, row.names=TRUE, sep="\t", file = "markers_out.tsv")
}

#+ echo = FALSE
if (showcode == TRUE && comparison == TRUE) print("Compare")
#+ echo = `showcode`, warning = `warn`, include = `comparison`
  Seurat::DefaultAssay(seuset) <- "ADT"
  seuset <- Seurat::NormalizeData(seuset, normalization.method = "CLR", margin = 2)
  Seurat::DefaultAssay(seuset) <- "RNA"
  seuset <- Seurat::NormalizeData(seuset, normalization.method = "CLR", margin = 2, assay = "ADT")
if (comparison == TRUE) {
  for(x in feature) {
    Seurat::DefaultAssay(seuset) <- "ADT"
    p1 <- Seurat::FeaturePlot(seuset, x, cols = c("lightgrey", "red")) + ggplot2::ggtitle(paste0("Protein:", " ", x))
    Seurat::DefaultAssay(seuset) <- "RNA"
    p2 <- Seurat::FeaturePlot(seuset, x) + ggplot2::ggtitle(paste0("RNA:", " ", x))
    print(p1 | p2)
    label <- as.character(paste0(Seurat::Key(seuset[["ADT"]]), x))
    print(Seurat::VlnPlot(seuset, paste0("rna_", x)))
    print(Seurat::VlnPlot(seuset, paste0("adt_", x)))
  }
}

#+ echo = FALSE
if (showcode == TRUE) print("Cite-seq")
#+ echo = `showcode`, warning = `warn`, include = `marker_compare`
rna_markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold, assay="RNA")
protein_markers <- Seurat::FindAllMarkers(seuset, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold, assay="ADT")
if (marker_compare == TRUE) {
  data.table::fwrite(x = rna_markers, sep="\t", file = "rna_out.tsv")
  data.table::fwrite(x = protein_markers, sep="\t", file = "protein_out.tsv")
}
toprna <- dplyr::top_n(dplyr::group_by(rna_markers, cluster), n=5, avg_log2FC)
toprna <- head(as.list(unique(as.data.frame(toprna)$gene)), top_x)
topprot <- dplyr::top_n(dplyr::group_by(protein_markers, cluster), n=5, avg_log2FC)
topprot <- head(as.list(unique(as.data.frame(topprot)$gene)), top_x)
if(marker_compare == TRUE) {
  pdf(file="citeseq_out.pdf")
  rna_labels <- as.vector(toprna)
  rna_labels <- rna_labels[!duplicated(rna_labels)]
  prot_labels <- as.vector(topprot)
  prot_labels <- prot_labels[!duplicated(prot_labels)]
  for(rnamarker in rna_labels) {
    rnamarker <-  paste("rna_", rnamarker, sep = "")
    for(protmarker in prot_labels) {
      protmarker <- paste("adt_", protmarker, sep="")
      plot <- Seurat::FeatureScatter(seuset, feature1 = rnamarker, feature2 = protmarker) + ggplot2::ggtitle(paste0(rnamarker, " vs ", protmarker))
      print(plot)
    }
  }
  for(rnamarker in rna_labels) {
    rnamarker <-  paste("rna_", rnamarker, sep = "")
    for(rnamarker2 in rna_labels) {
      rnamarker2 <- paste("rna_", rnamarker2, sep="")
      plot <- Seurat::FeatureScatter(seuset, feature1 = rnamarker, feature2 = rnamarker2) + ggplot2::ggtitle(paste0(rnamarker, " vs ", rnamarker2))
      print(plot)
    }
  }
  for(protmarker in prot_labels) {
    protmarker <-  paste("adt_", protmarker, sep = "")
    for(protmarker2 in prot_labels) {
      protmarker2 <- paste("adt_", protmarker2, sep="")
      plot <- Seurat::FeatureScatter(seuset, feature1 = protmarker, feature2 = protmarker2) + ggplot2::ggtitle(paste0(protmarker, " vs ", protmarker2))
      print(plot)
    }
  }
  dev.off()
}

# nolint end
