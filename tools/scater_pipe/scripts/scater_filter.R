#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

# Load libs, common functions, source RaceID and Galaxy Params
source(paste(script_dir, "common.R", sep="/"))

# Read input data
x <- read.table(count_matrix, header=col_header, comment=commchar, row.names=row_head)

#message("Count matrix with %.0f cells and %.0f genes", dim(x)[1], dim(x)[2])

# Perform the spike-in filtering
x_filt <- x
if (filterout != "" ){
    x_filt <- x[!grepl(filterout,x_filt),]
}

sce <- SingleCellExperiment(assays = list(counts = as.matrix(x_filt)))
sce <- calculateQCMetrics(sce)
cnts <- counts(sce)

sce_f <- sce[rowSums(cnts) > min_gene_trans, colSums(cnts) > min_cell_trans]

plotter("OUT_prefilter_totalcounts", function(){
    hist(
        sce$total_counts, breaks=100,
        xlab="Total Transcripts per Cell",
        ylab="Frequency",
        main=""
    )
    abline(v=min_gene_trans, col='red')
})

plotter("OUT_prefilter_totalgenes" , function(){
    hist(sce$total_features, breaks=100,
         xlab="Total Gene Transcripts per Cell",
         ylab="Frequency",
         main="")
    abline(v=min_gene_trans, col = 'red')
})


# Table of genes with the number of cells having a transcript
# count for that gene higher than gdetect
cell_counts_per_gene <- rowSums(counts(sce_f) > gdetect_mintrans)
list_of_genes_with_n_cells <- cell_counts_per_gene > gdetect_mincells

gdetect_matrix <- sce_f[list_of_genes_with_n_cells,]

tab <- rbind(
    dim(x),
    dim(sce),
    dim(sce_f),
    dim(gdetect_matrix)
)
colnames(tab) <- c("Genes", "Cells")
rownames(tab) <- c("Initial Matrix", "Filtered names", "Filtered counts", "Filtered Gene Detection")

write.table(tab, "table.stats")
write.table(counts(gdetect_matrix), "table.counts")
saveRDS(gdetect_matrix, "table.rds")
