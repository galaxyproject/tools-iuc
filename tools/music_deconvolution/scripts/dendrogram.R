
suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to generate a clustering dendrogram of cell types, using the prior labelling from scRNA.

bulk_eset = readRDS('https://xuranw.github.io/MuSiC/data/Mousebulkeset.rds')
scrna_eset = readRDS('https://xuranw.github.io/MuSiC/data/Mousesubeset.rds')

print(bulk_eset)
print(scrna_eset)
# Download EMTAB single cell dataset from Github

celltypes = c('Endo', 'Podo', 'PT', 'LOH', 'DCT', 'CD-PC', 'CD-IC',
              'Fib', 'Macro', 'Neutro','B lymph', 'T lymph', 'NK')
clusters_label = 'cellType'
samples_label = 'sampleID'
list.of.markers = read.table() ## load('../IEmarkers.RData')
## how can we translate this to IEmarkers?

clusters.type = list(C1 = 'Neutro', C2 = 'Podo',
                     C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'),
                     C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))
## Maybe a table with the following groups:
## C1 Neutro 
## C2 Podo
## C3 Endo CD-PC LOH CD-IC DCT PT
## C4 Macro Fib  B lymph   NK   T lymph

## the idea is that the user explores the dataset to 
levels(scrna_eset$cellType)
## [1] "Endo"     "Podo"     "PT"       "LOH"      "DCT"      "CD-PC"    "CD-IC"    "CD-Trans" "Novel1"  
##[10] "Fib"      "Macro"    "Neutro"   "B lymph"  "T lymph"  "NK"       "Novel2"

## Perform the estimation
## Produce the first step information
sub.basis = music_basis(scrna_eset, clusters = clusters_label,
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


## We then perform bulk tissue cell type estimation with pre-grouping of cell types.
clusters.type = list(C1 = 'Neutro', C2 = 'Podo',
                     C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'),
                     C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))

cl.type = as.character(scrna_eset$cellType)

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(scrna_eset)$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse
#       C1        C2       C31       C32       C33       C34       C35       C36       C41       C42 
# "Neutro"    "Podo"    "Endo"   "CD-PC"     "LOH"   "CD-IC"     "DCT"      "PT"   "Macro"     "Fib" 
#      C43       C44       C45 
##"B lymph"      "NK" "T lymph"


## We load cell markers
# This RData file provides two vectors of gene names Epith.marker and Immune.marker

# We now construct the list of group marker
IEmarkers = list(NULL, NULL, Epith.marker, Immune.marker)
names(IEmarkers) = c('C1', 'C2', 'C3', 'C4')
# The name of group markers should be the same as the cluster names

Est.mouse.bulk = music_prop.cluster(bulk.eset = bulk_eset, sc.eset = scrna_eset,
                                    group.markers = IEmarkers, clusters = 'cellType',
                                    groups = 'clusterType', samples = 'sampleID',
                                    clusters.type = clusters.type)

