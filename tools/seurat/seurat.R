options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    make_option(c("-counts","--counts"), type="character", help="Counts file"),
    make_option(c("-numPCs","--numPCs"), type="integer", help="Number of PCs to use in plots"),
    make_option(c("-min.cells","--min.cells"), type="integer", help="Minimum cells to include"),
    make_option(c("-min.genes","--min.genes"), type="integer", help="Minimum genes to include"),
    make_option(c("-low.thresholds","--low.thresholds"), type="double", help="Low threshold for filtering cells"),
    make_option(c("-high.thresholds","--high.thresholds"), type="double", help="High threshold for filtering cells"),
    make_option(c("-x.low.cutoff","--x.low.cutoff"), type="double", help="X-axis low cutoff for variable genes"),
    make_option(c("-x.high.cutoff","--x.high.cutoff"), type="double", help="X-axis high cutoff for variable genes"),
    make_option(c("-y.cutoff","--y.cutoff"), type="double", help="Y-axis cutoff for variable genes"),
    make_option(c("-cells.use","--cells.use"), type="integer", help="Cells to use for PCHeatmap"),
    make_option(c("-resolution","--resolution"), type="double", help="Resolution in FindClusters"),
    make_option(c("-min.pct","--min.pct"), type="double", help="Minimum percent cells in FindClusters"),
    make_option(c("-logfc.threshold","--logfc.threshold"), type="double", help="LogFC threshold in FindClusters"),
    make_option(c("-rds","--rds"), type="logical", help="Output Seurat RDS object")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

counts <- read.delim(args$counts, row.names=1)
seuset <- CreateSeuratObject(raw.data = counts, min.cells = args$min.cells, min.genes = args$min.cells)

# Open PDF for plots
pdf("out.pdf")

VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 2)
GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")

print("Filtering cells")
if (!is.null(args$low.thresholds)){
    lowthresh <- args$low.thresholds
} else {
    lowthresh <- "-Inf"
}
if (!is.null(args$high.thresholds)){
    highthresh <- args$high.thresholds
} else {
    highthresh <- "Inf"
}
seuset <- FilterCells(object = seuset, subset.names = c("nUMI"), 
    low.thresholds=c(lowthresh), high.thresholds = c(highthresh))

print("Normalizing the data")
seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", 
    scale.factor = 10000)

print("Finding variable genes")
seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = args$x.low.cutoff, 
    x.high.cutoff = args$x.high.cutoff,,
    y.cutoff = args$y.cutoff
)

print("Scaling the data and removing unwanted sources of variation")
seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI"))

print("Performing PCA analysis")
seuset <- RunPCA(object = seuset, pc.genes = seuset@var.genes)
VizPCA(object = seuset, pcs.use = 1:2)
PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
PCHeatmap(
    object = seuset, 
    pc.use = 1:args$numPCs, 
    cells.use = args$cell.use, 
    do.balanced = TRUE, 
    label.columns = FALSE,
    use.full = FALSE
)

print("Determining statistically significant principal components")
seuset <- JackStraw(object = seuset, num.replicate = 100, display.progress= FALSE)
JackStrawPlot(object = seuset, PCs = 1:args$numPCs)
PCElbowPlot(object = seuset)

print("Clustering the cells")
seuset <- FindClusters(
    object = seuset, 
    reduction.type = "pca", 
    dims.use = 1:args$numPCs, 
    resolution = args$resolution,
    print.output = 0, 
    save.SNN = TRUE
)

print("Running non-linear dimensional reduction (tSNE)")
seuset <- RunTSNE(object = seuset, dims.use = 1:args$numPCs, do.fast = TRUE)
TSNEPlot(object = seuset)

print("Finding differentially expressed genes (cluster biomarkers)")
markers <- FindAllMarkers(object = seuset, only.pos = TRUE, min.pct = args$min.pct,
    logfc.threshold = args$logfc.threshold)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seuset, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# Close PDF for plots
dev.off()

if (!is.null(args$rds) ) {
  saveRDS(seuset, "Seurat.rds")
}

sessionInfo()