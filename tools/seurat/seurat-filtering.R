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
    make_option("--counts", type="character", help="Counts file"),
    make_option("--min.cells", type="integer", help="Minimum cells that detected the gene in ordrer to be included"),
    make_option("--min.genes", type="integer", help="Minimum genes detected by the cell to be included"),
    make_option("--low.thresholds", default=NULL, type="double", help="Low threshold for filtering cells"),
    make_option("--high.thresholds", default=NULL, type="double", help="High threshold for filtering cells"),
    make_option("--rds", type="character", help="Output Seurat RDS object"),
    make_option("--pdf", type="character", help="Output PDF file for plots")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

counts <- read.delim(args$counts, row.names=1)
seuset <- CreateSeuratObject(raw.data = counts, min.cells = args$min.cells, min.genes = args$min.genes)

# Open PDF for plots
pdf(args$pdf)

VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 2)
GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene", col=rgb(0,0,0, alpha=0.5))

if (!is.null(args$low.thresholds)){
    lowthresh <- args$low.thresholds
    abline(v=lowthresh, col="red", lwd=2, lty=2)
} else {
    lowthresh <- -Inf
}

if (!is.null(args$high.thresholds)){
    highthresh <- args$high.thresholds
    abline(v=highthresh, col="red", lwd=2, lty=2)
} else {
    highthresh <- Inf
}
seuset <- FilterCells(object = seuset, subset.names = c("nUMI"), 
    low.thresholds=c(lowthresh), high.thresholds = c(highthresh))

# Close PDF for plots
dev.off()

saveRDS(seuset, args$rds)

sessionInfo()
