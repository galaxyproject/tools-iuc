##
suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to generate a clustering dendrogram of cell types, using the prior labelling from scRNA.

read_list <- function(lfile){
    if (lfile=="None"){
        return(NULL)
    }
    return(read.table(file=lfile, header=F, stringsAsFactors=F)$V1)
}

args = commandArgs(trailingOnly=TRUE)
source(args[1])

## bulk_eset = readRDS('../test-data/Mousebulkeset.rds')
## scrna_eset = readRDS('../test-data/Mousesubeset.rds')

## ## The first user inspects the scrna dataset to find these types
## celltypes = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC', 'Fib', 'Macro', 'Neutro','B lymph', 'T lymph', 'NK')
## ## "Endo,Podo,PT,LOH,DCT,CD-PC,CD-IC,Fib,Macro,Neutro,B lymph,T lymph,NK"
## celltypes_label = 'cellType'
## clustertype_label = 'clusterType'
## samples_label = 'sampleID'
## ## We then perform bulk tissue cell type estimation with pre-grouping of cell types.
## ## CN, list_of_cell_types, marker genes name, marker genes list
## ## e.g.
## ##   C1, c("Neutro"), NULL, NULL
## ##   C2, c("Podo"), NULL, NULL
## ##   C3, c("Endo","CD-PC","LOH","CD-IC","DCT","PT"), Epith.marker, <data>
## ##   C4, c("Macro","Fib","B lymph","NK","T lymph"), Immune.marker, <data>

## data.to.use = list(
##     "C1" = list(cell.types = c("Neutro"),
##                 marker.names=NULL,
##                 marker.list=NULL),
##     "C2" = list(cell.types = c("Podo"),
##                 marker.names=NULL,
##                 marker.list=NULL),
##     "C3" = list(cell.types = c("Endo","CD-PC","LOH","CD-IC","DCT","PT"),
##                 marker.names = "Epithelial",
##                 marker.list = read_list("../test-data/epith.markers")),
##     "C4" = list(cell.types = c("Macro","Fib","B lymph","NK","T lymph"),
##                 marker.names = "Immune",
##                 marker.list = read_list("../test-data/immune.markers"))
## )
grouped_celltypes = lapply(data.to.use, function(x){ x$cell.types })
marker_groups = lapply(data.to.use, function(x){ x$marker.list })
names(marker_groups) = names(data.to.use) 


## Perform the estimation
## Produce the first step information
sub.basis = music_basis(scrna_eset, clusters = celltypes_label,
                        samples = samples_label,
                        select.ct = celltypes)

## Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(sub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
## Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
## Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(sub.basis$M.theta + 1e-8)), method = "euclidean")
## Hierarchical clustering using Complete Linkage
hc2 <- hclust(d, method = "complete")
## Plot the obtained dendrogram
pdf(file=outfile_pdf, width=8, height=8)
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

cl.type = as.character(scrna_eset[[celltypes_label]])

for(cl in 1:length(grouped_celltypes)){
  cl.type[cl.type %in% grouped_celltypes[[cl]]] = names(grouped_celltypes)[cl]
}
pData(scrna_eset)[[clustertype_label]] = factor(cl.type, levels = c(names(grouped_celltypes), 'CD-Trans', 'Novel1', 'Novel2'))

Est.bulk = music_prop.cluster(bulk.eset = bulk_eset, sc.eset = scrna_eset,
                              group.markers = marker_groups, clusters = celltypes_label,
                              groups = clustertype_label, samples = samples_label,
                              clusters.type = grouped_celltypes)

print(Est.bulk)

dev.off()