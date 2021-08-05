## 
suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to generate a clustering dendrogram of cell types, using the prior labelling from scRNA.

bulk_eset = readRDS('../test-data/Mousebulkeset.rds')
scrna_eset = readRDS('../test-data/Mousesubeset.rds')

print(bulk_eset)
print(scrna_eset)
# Download EMTAB single cell dataset from Github

## The first user inspects the scrna dataset to find these types
celltypes = c('Endo', 'Podo', 'PT', 'LOH', 'DCT',
              'CD-PC', 'CD-IC', 'Fib', 'Macro',
              'Neutro','B lymph', 'T lymph', 'NK')
cell_type_label = 'cellType'
cluster_type_label = 'clusterType'
samples_label = 'sampleID'
## We then perform bulk tissue cell type estimation with pre-grouping of cell types.
tab <- read.table('../test-data/cluster_groups.tsv', sep="\t", header=F, row.names=1, stringsAsFactors=FALSE)
pre_grouped_cell_types = lapply(split(tab, row.names(tab)),
                       function(x){
                           unlist(strsplit(x$V2, split=","))
                       })
## Perform the estimation
## Produce the first step information
sub.basis = music_basis(scrna_eset, clusters = cell_type_label,
                        samples = samples_label, 
                        select.ct = celltypes)

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(sub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(sub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

cl.type = as.character(scrna_eset[[cell_type_label]])

for(cl in 1:length(pre_grouped_cell_types)){
  cl.type[cl.type %in% pre_grouped_cell_types[[cl]]] = names(pre_grouped_cell_types)[cl]
}
pData(scrna_eset)[[cluster_type_label]] = factor(cl.type, levels = c(names(pre_grouped_cell_types), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
##cell_types_flattened = unlist(pre_grouped_cell_types)
##cell_types_flattened
#       C1        C2       C31       C32       C33       C34       C35       C36       C41       C42 
# "Neutro"    "Podo"    "Endo"   "CD-PC"     "LOH"   "CD-IC"     "DCT"      "PT"   "Macro"     "Fib" 
#      C43       C44       C45 
##"B lymph"      "NK" "T lymph"



## We load cell markers
## This RData file provides two vectors of gene names Epith.marker and Immune.marker
load('../test-data/IEmarkers.RData')
##


# We now construct the list of group marker
IEmarkers = list(NULL, NULL, Epith.marker, Immune.marker)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4')
# The name of group markers should be the same as the cluster names

Est.mouse.bulk = music_prop.cluster(bulk.eset = bulk_eset, sc.eset = scrna_eset,
                                    group.markers = IEmarkers, clusters = cell_type_label,
                                    groups = cluster_type_label, samples = samples_label,
                                    clusters.type = pre_grouped_cell_types)

