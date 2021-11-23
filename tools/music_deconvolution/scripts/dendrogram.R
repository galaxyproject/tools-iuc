##
suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
## We use this script to generate a clustering dendrogram of cell
## types, using the prior labelling from scRNA.

read_list <- function(lfile) {
    if (lfile == "None") {
        return(NULL)
    }
    return(read.table(file = lfile, header = FALSE,
                      stringsAsFactors = FALSE)$V1)
}

args <- commandArgs(trailingOnly = TRUE)
source(args[1])


## Perform the estimation
## Produce the first step information
sub.basis <- music_basis(scrna_eset, clusters = celltypes_label,
                         samples = samples_label,
                         select.ct = celltypes)

## Plot the dendrogram of design matrix and cross-subject mean of
## realtive abundance
## Hierarchical clustering using Complete Linkage
d1 <- dist(t(log(sub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
hc1 <- hclust(d1, method = "complete")
## Hierarchical clustering using Complete Linkage
d2 <- dist(t(log(sub.basis$M.theta + 1e-8)), method = "euclidean")
hc2 <- hclust(d2, method = "complete")


if (length(data.to.use) > 0) {
    ## We then perform bulk tissue cell type estimation with pre-grouping
    ## of cell types: C, list_of_cell_types, marker genes name, marker
    ## genes list.
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
    grouped_celltypes <- lapply(data.to.use, function(x) {
        x$cell.types
    })
    marker_groups <- lapply(data.to.use, function(x) {
        x$marker.list
    })
    names(marker_groups) <- names(data.to.use)


    cl_type <- as.character(scrna_eset[[celltypes_label]])

    for (cl in seq_len(length(grouped_celltypes))) {
        cl_type[cl_type %in%
                grouped_celltypes[[cl]]] <- names(grouped_celltypes)[cl]
    }
    pData(scrna_eset)[[clustertype_label]] <- factor(
        cl_type, levels = c(names(grouped_celltypes),
                            "CD-Trans", "Novel1", "Novel2"))

    est_bulk <- music_prop.cluster(
        bulk.eset = bulk_eset, sc.eset = scrna_eset,
        group.markers = marker_groups, clusters = celltypes_label,
        groups = clustertype_label, samples = samples_label,
        clusters.type = grouped_celltypes
    )

    estimated_music_props <- est_bulk$Est.prop.weighted.cluster
    ## NNLS is not calculated here

    ## Show different in estimation methods
    ## Jitter plot of estimated cell type proportions
    methods_list <- c("MuSiC")

    jitter_fig <- Jitter_Est(
        list(data.matrix(estimated_music_props)),
        method.name = methods_list, title = "Jitter plot of Est Proportions",
        size = 2, alpha = 0.7) +
        theme_minimal() +
        labs(x = element_blank(), y = element_blank()) +
        theme(axis.text = element_text(size = 6),
              axis.text.x = element_blank(),
              legend.position = "none")

    plot_box <- Boxplot_Est(list(
        data.matrix(estimated_music_props)),
        method.name = methods_list) +
        theme_minimal() +
        labs(x = element_blank(), y = element_blank()) +
        theme(axis.text = element_text(size = 6),
              axis.text.x = element_blank(),
              legend.position = "none")

    plot_hmap <- Prop_heat_Est(list(
        data.matrix(estimated_music_props)),
        method.name = methods_list) +
        labs(x = element_blank(), y = element_blank()) +
        theme(axis.text.y = element_text(size = 6),
              axis.text.x = element_text(angle = -90, size = 5),
              plot.title = element_text(size = 9),
              legend.key.width = unit(0.15, "cm"),
              legend.text = element_text(size = 5),
              legend.title = element_text(size = 5))

}
    
pdf(file = outfile_pdf, width = 8, height = 8)
par(mfrow = c(1, 2))
plot(hc1, cex = 0.6, hang = -1, main = "Cluster log(Design Matrix)")
plot(hc2, cex = 0.6, hang = -1, main = "Cluster log(Mean of RA)")
if (length(data.to.use) > 0) {
    plot_grid(jitter_fig, plot_box, plot_hmap, ncol = 2, nrow = 2)
}
message(dev.off())

if (length(data.to.use) > 0) {
    write.table(estimated_music_props,
                file = outfile_tab, quote = F, col.names = NA, sep = "\t")
}
