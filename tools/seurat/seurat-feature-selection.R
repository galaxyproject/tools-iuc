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
    make_option("--data", type="character", help="Counts file"),
    make_option("--x.low.cutoff", type="double", default=0.1, help="Low threshold for average expression"),
    make_option("--x.high.cutoff", type="double", default=8, help="High threshold for average expression"),
    make_option("--y.cutoff", type="double", default=1, help="Dispersion cutoff"),
    make_option("--selection.method", type="character", help="How to select high variable genes ? Based on dispersion value or with the different cutoff"),
    make_option("--top.genes", type="numeric", help="Number of genes to keep"),
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

seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = args$x.low.cutoff, 
    x.high.cutoff = args$x.high.cutoff,,
    y.cutoff = args$y.cutoff,
    selection.method = args$selection.method,
    top.genes = args$top.genes
)

# Close PDF for plots
dev.off()

save(seuset, file = args$rds)

print("Number of feature selected: ")
print(length(seuset@var.genes))

sessionInfo()
