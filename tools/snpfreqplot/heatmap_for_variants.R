#!/usr/bin/env R

library(pheatmap)
library(RColorBrewer)
library(tidyverse)


                                        # Galaxy Variables
samples = data.frame(ids=c("436", "437", "438", "439", "440", "441", "442", "443", "444"),
                     ## the full file paths below are populated by the galaxy wrapper
                     files=c("./input//Galaxy436-[SnpSift_Extract_Fields_on_data_367].tabular",
                             "./input//Galaxy437-[SnpSift_Extract_Fields_on_data_370].tabular",
                             "./input//Galaxy438-[SnpSift_Extract_Fields_on_data_373].tabular",
                             "./input//Galaxy439-[SnpSift_Extract_Fields_on_data_376].tabular",
                             "./input//Galaxy440-[SnpSift_Extract_Fields_on_data_379].tabular",
                             "./input//Galaxy441-[SnpSift_Extract_Fields_on_data_382].tabular",
                             "./input//Galaxy442-[SnpSift_Extract_Fields_on_data_385].tabular",
                             "./input//Galaxy443-[SnpSift_Extract_Fields_on_data_388].tabular",
                             "./input//Galaxy444-[SnpSift_Extract_Fields_on_data_391].tabular"))
##
date.has <- FALSE                   ## do the file names contain a date (format: dd.mm.yyyy) and you want to sort
date.regex <- "\\d{2}[[:punct:]]\\d{2}[[:punct:]]\\d{4}"
variant.frequency <- 0.1                   ## adjust the variant frequency
img.multiplier_width <- 0.5                ## pdf weight multiplier - adjust if needed
img.multiplier_height <- 0.6               ## pdf height multiplier - adjust if needed
brewer.color_gene_annotation <- "Set3"     ## adjust color of gene annotation ("Set2" or "Paired")
pheat.clustering <- FALSE                  ## should the samples be clustered?
pheat.clustering_method <- "ward.D2"       ## what clustering method -> see ?hclust for further information
pheat.number_of_clusters <- 5              ## do you assume a particular amount of clusters?
out.filepdf <- "output/Heatmap.pdf"



if (date.has){
    samples$dates <- as.Date(regmatches(samples$files,
                                        regexpr(date.regex, samples$files)),
                             format="%d.%m.%y")
    samples$ids <- factor(samples$ids,
                          ## files are processed in factor order if date is given
                          levels=samples$ids[order(sort(samples$dates))])
}

extractAll_data <- function(id){
    file <- (samples %>% filter(ids==id))$files
    variants <- read.table(file, header = T, sep = "\t")
    tmp <- variants %>% mutate(posalt=paste(POS,ALT)) %>% select(posalt, AF)
    colnames(tmp) <- c("Mutation", id)
    return(tmp)
}

extractAll_annots <- function(id){
    file <- (samples %>% filter(ids==id))$files
    variants <- read.table(file, header = T, sep = "\t")
    tmp <- variants %>% mutate(posalt=paste(POS,ALT),
                               effect=EFF....EFFECT,
                               gene=EFF....GENE) %>%
        select(posalt, effect, gene)
}


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
pdf(out.filepdf,
    width = img.multiplier_width*ncol(final),      ## the pdf scales with the number
    height = img.multiplier_height*nrow(final))    ## of variants and samples

pheatmap(final,
         color = my_colors(100),
         cellwidth = 20,
         cellheight = 20,
         clustering_method = pheat.clustering_method,
         cluster_rows = pheat.clustering,
         cluster_cols = F,
         cutree_rows = pheat.number_of_clusters,
         annotation_col = ann_final,
         annotation_colors = color_list,
         gaps_col = gap_vector)

dev.off()
