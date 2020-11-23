#!/usr/bin/env R

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))

extractAll_data <- function(id){
    variants = readAndProcess(id)
    tmp = variants %>% mutate(posalt=uni.select) %>%
        select(posalt, AF)
    colnames(tmp) <- c("Mutation", id)
    return(tmp)
}

extractAll_annots <- function(id){
    variants = readAndProcess(id)
    tmp <- variants %>%
        mutate(posalt=paste(POS,ALT),
               effect=EFF....EFFECT, gene=EFF....GENE) %>%
        select(posalt, effect, gene)
    return(tmp)
}
                                        # M A I N
stopifnot(exists("samples"))
                                        # process allele frequencies
processed.files <- lapply(samples$ids, extractAll_data)
final <- as_tibble(
    processed.files %>%
    reduce(full_join, by = "Mutation", copy=T))

final <- final[str_order(final$Mutation, numeric = T),] %>%
    column_to_rownames("Mutation")                          ## sort and set rownames
final[final < variant.frequency] <- NA                      ## adjust the variant frequency:
final <- final[rowSums(is.na(final)) !=ncol(final), ]
final <- t(final)
final[is.na(final)] <- 0
class(final) <- "numeric"

                                        # add annotations
processed.annots <-lapply(samples$ids, extractAll_annots)   ## readout annotations
ann_final <- processed.annots %>%
    reduce(function(x,y){unique(rbind(x,y))}) %>%
    filter(posalt %in% colnames(final))                     ## apply frequency filter
ann_final <- as_tibble(ann_final[str_order(ann_final$posalt, numeric = T),]) %>%
    column_to_rownames("posalt")                            ## sort

                                        # rename annotations
ann_final$gene <- sub("^$", "NCR", ann_final$gene)
ann_final$effect <- sub("^$", "non-coding", ann_final$effect)
ann_final$effect[ann_final$effect=="NON_SYNONYMOUS_CODING"] <-"non-syn"
ann_final$effect[ann_final$effect=="SYNONYMOUS_CODING"] <- "syn"
ann_final$effect[ann_final$effect=="CODON_CHANGE_PLUS_CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="FRAME_SHIFT"] <- "frame shift"
ann_final$effect[ann_final$effect=="STOP_GAINED"] <- "stop gained"
ann_final$effect[ann_final$effect=="CODON_DELETION"] <- "deletion"
ann_final$effect[ann_final$effect=="FRAME_SHIFT+STOP_GAINED"] <- "stop gained"
ann_final$effect[ann_final$effect=="CODON_CHANGE_PLUS_CODON_INSERTION"] <- "insertion"
ann_final$effect[ann_final$effect=="INSERTION"] <- "insertion"

                                        # automatically determine gaps for the heatmap
gap_vector <- which(!(ann_final$gene[1:length(ann_final$gene)-1] ==
                      ann_final$gene[2:length(ann_final$gene)]))

                                        # colormanagement
my_colors <- colorRampPalette(c("grey93","brown","black"))  #colormangment heatmap
count <- length(unique(ann_final$gene))                     #colormanagement annotations (genes)
gene_color <- c(brewer.pal(brewer.color_gene_annotation, n=count))
names(gene_color) = unique(ann_final$gene)

                                        # colormanagement annotations (effect)
colors <- c()
color_test <- function(eff, col, colors){
    if(eff %in% ann_final$effect) colors <<- c(colors,col)
}
color_test("non-coding", "white", colors)
color_test("syn", "green", colors)
color_test("non-syn", "orange", colors)
color_test("deletion", "red", colors)
color_test("frame shift", "black", colors)
color_test("stop gained", "grey", colors)
color_test("insertion", "blue", colors)
##
all_colors <- data.frame(
    color=c("white", "green", "orange", "red", "black", "grey", "blue"),
    name = c("non-coding", "syn", "non-syn", "deletion", "frame shift", "stop gained", "insertion"))
subset_colors <- subset(all_colors, color %in% colors)
effect_color <- subset_colors$color
names(effect_color) = subset_colors$name
color_list <- list(gene_color = gene_color, effect_color = effect_color)
names(color_list) <- c("gene", "effect")

                                        # visualize heatmap
 if (pheat.number_of_clusters > length(samples$ids)){
    print(paste0("[INFO] Number of clusters: User-specified clusters (",
                 pheat.number_of_clusters,
                 ") is greater than the number of samples (",
                 length(samples$ids), ")"))
    pheat.number_of_clusters = length(samples$ids)
    print(paste0("[INFO] Number of clusters: now set to ",
                 pheat.number_of_clusters))
}

## image trickery needed here:
## - the heatmap file type is determined by the extension, but
##   galaxy only serves .dat files, so we need to create a
##   a new output file with that extension, and move to the
##   the galaxy output dat file.

out.tmpfile = tempfile(fileext=paste0(".", tolower(out.type)))

pheatmap(final,
         color = my_colors(100),
         cellwidth = cell.width,
         cellheight = cell.height,
         clustering_method = pheat.clustering_method,
         cluster_rows = pheat.clustering,
         cluster_cols = F,
         cutree_rows = pheat.number_of_clusters,
         annotation_col = ann_final,
         annotation_colors = color_list,
         filename = out.tmpfile,
         gaps_col = gap_vector)

file.rename(from=out.tmpfile, to=out.file)
