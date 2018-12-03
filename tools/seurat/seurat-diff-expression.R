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
    make_option("--data", type="character", help="Seurat RDS object"),
    make_option("--min.pct", type="double", help="Minimum percent cells in FindAllMarkers"),
    make_option("--logfc.threshold", type="double", help="LogFC threshold in FindAllMarkers"),
    make_option("--rds", type="character", help="Output Seurat RDS object"),
    make_option("--tab", type="character", help="Output table of differentially expressed genes"),
    make_option("--pdf", type="character", help="Output PDF file for plots")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

load(args$data)


#Check if there is the PCA slot
if(!("pca" %in% names(seuset@dr))) stop("You need to perform a PCA before the clustering.")

# Open PDF for plots
pdf(args$pdf)


print("Finding differentially expressed genes (cluster biomarkers)")
markers <- FindAllMarkers(object = seuset, only.pos = TRUE, min.pct = args$min.pct,
    logfc.threshold = args$logfc.threshold)
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = seuset, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

topgenes <- markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

#check if tsne available
if(("tsne" %in% names(seuset@dr))){
    dim.red = "tsne"
} else {
    dim.red = "pca"
}

FeaturePlot(object = seuset, features.plot = topgenes$gene, cols.use = c("grey", "blue"), 
    reduction.use = dim.red)

# Close PDF for plots
dev.off()

save(seuset, file = args$rds)
write.table(markers, args$tab, sep="\t", quote=F, row.names = F)

sessionInfo()
