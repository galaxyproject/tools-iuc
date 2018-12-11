options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dplyr)
    library(gdata)
    library(optparse)
})

option_list <- list(
    make_option("--data", type="character", help="Seurat RDS object"),
    make_option("--numPCs", type="integer", help="Number of PCs to use in plots"),
    make_option("--cells.use", type="integer", default=NULL, help="Cells to use for PCHeatmap"),
    make_option("--col.pca", type="character", default="orig.ident", help="Variable used for PCA color"),
    make_option("--do.tsne", type="logical", help="Proceed t-SNE"),
    make_option("--numPCs.tsne", type="integer", help="Number of PCs to use in plots for t-SNE"),
    make_option("--col.tsne", type="character", default="orig.ident", help="Variable used for t-SNE color"),
    make_option("--rds", type="character", help="Output Seurat RDS object"),
    make_option("--pdf", type="character", help="Output PDF file for plots")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

a <- load(args$data)
#change imported object name in seuset if it's not the case
if(!exists("seuset")) mv(a, "seuset")

# Open PDF for plots
pdf(args$pdf)

print("Performing PCA analysis")
seuset <- RunPCA(object = seuset, pc.genes = seuset@var.genes)
VizPCA(object = seuset, pcs.use = 1:args$numPCs)
PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2, group.by=args$col.pca)
PCHeatmap(
    object = seuset, 
    pc.use = 1:args$numPCs, 
    cells.use = args$cells.use, 
    do.balanced = TRUE, 
    label.columns = FALSE,
    use.full = FALSE
)

print("Determining statistically significant principal components")
seuset <- JackStraw(object = seuset, num.replicate = 100, display.progress= FALSE)
JackStrawPlot(object = seuset, PCs = 1:args$numPCs)
PCElbowPlot(object = seuset)

if(args$do.tsne == T){
    print("Running non-linear dimensional reduction (tSNE)")
    seuset <- RunTSNE(object = seuset, dims.use = 1:args$numPCs.tsne, do.fast = TRUE)
    TSNEPlot(object = seuset, group.by=args$col.tsne)
}

# Close PDF for plots
dev.off()

save(seuset, file = args$rds)


sessionInfo()
