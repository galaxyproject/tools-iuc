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
    make_option("--numPCs", type="integer", help="Number of PCs to use in FindClusters"),
    make_option("--resolution", type="double", help="Resolution in FindClusters"),
    make_option("--rds", type="character", help="Output Seurat RDS object"),
    make_option("--pdf", type="character", help="Output PDF file for plots")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

a <- load(args$data)
#change imported object name in seuset if it's not the case
if(!exists("seuset")) mv(a, "seuset")

#Check if there is the PCA slot
if(!("pca" %in% names(seuset@dr))) stop("You need to perform a PCA before the clustering.")

# Open PDF for plots
pdf(args$pdf)

print("Clustering the cells")
seuset <- FindClusters(
    object = seuset, 
    dims.use = 1:args$numPCs, 
    resolution = args$resolution,
    print.output = 0, 
    save.SNN = TRUE
)

PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
if(!("tsne" %in% names(seuset@dr))){
    TSNEPlot(object = seuset)
}

# Close PDF for plots
dev.off()

save(seuset, file = args$rds)


sessionInfo()
