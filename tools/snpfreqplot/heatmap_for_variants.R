#!/usr/bin/env R

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))

fapply <- function(vect_ids, func) {
    #' List apply but preserve the names
    res <- lapply(vect_ids, func)
    names(res) <- vect_ids
    return(res)
}

                                        # M A I N
stopifnot(exists("samples"))
variant_files <- fapply(samples$ids, read_and_process)  # nolint

extractall_data <- function(id) {
    variants <- variant_files[[id]]
    tmp <- variants %>%
        mutate(unique_selectors = group_select) %>%
        select(unique_selectors, AF)
    colnames(tmp) <- c("Mutation", id)
    return(tmp)
}

extractall_annots <- function(id) {
    variants <- variant_files[[id]]
    tmp <- variants %>%
        mutate(unique_selectors = group_select,
               effect = EFF....EFFECT, gene = EFF....GENE) %>%
        select(unique_selectors, effect, gene)
    # allow "." as an alternative missing value in EFF.EFFECT and EFF.GENE
    tmp$effect <- sub("^\\.$", "", tmp$effect)
    tmp$gene <- sub("^\\.$", "", tmp$gene)
    return(tmp)
}

                                        # process allele frequencies
processed_files <- fapply(samples$ids, extractall_data)
final <- as_tibble(
    processed_files %>%
    reduce(full_join, by = "Mutation", copy = T))

final <- final[str_order(final$Mutation, numeric = T), ] %>%
    column_to_rownames("Mutation")              ## sort and set rownames
final[final < variant_frequency] <- NA          ## adjust the variant frequency:
final <- final[rowSums(is.na(final)) != ncol(final), ]
final <- t(final)
final[is.na(final)] <- 0
class(final) <- "numeric"

                                        # add annotations
## readout annotations
processed_annots <- fapply(samples$ids, extractall_annots)
ann_final <- processed_annots %>%
    reduce(function(x, y) {
        unique(rbind(x, y))}) %>%
    ## apply frequency filter
    filter(unique_selectors %in% colnames(final))
ann_final <- as_tibble(ann_final[str_order(
    ann_final$unique_selectors, numeric = T), ]) %>%
    column_to_rownames("unique_selectors")  ## sort

                                        # rename annotations
trans <- function(x, mapping, replace_missing=NULL) {
    # helper function for translating effects
    mapped <- mapping[[x]]
    if (is.null(mapped)) {
        if (is.null(replace_missing)) x else replace_missing
    } else {
        mapped
    }
}

# handle translation of classic SnpEff effects to sequence ontology terms
# The following list defines the complete mapping between classic and So effect
# terms even if not all of these are likely to appear in viral variant data.
classic_snpeff_effects_to_so <- list(
    "coding_sequence_variant", "coding_sequence_variant", "disruptive_inframe_deletion", "disruptive_inframe_insertion", "inframe_deletion", "inframe_insertion", "downstream_gene_variant", "exon_variant", "exon_loss_variant", "frameshift_variant", "gene_variant", "intergenic_variant", "intergenic_region", "conserved_intergenic_variant", "intragenic_variant", "intron_variant", "conserved_intron_variant", "missense_variant", "rare_amino_acid_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "5_prime_UTR_premature_start_codon_variant", "start_lost", "stop_gained", "stop_lost", "synonymous_variant", "start_retained_variant", "stop_retained_variant", "transcript_variant", "upstream_gene_variant", "3_prime_UTR_truncation_+_exon_loss_variant", "3_prime_UTR_variant", "5_prime_UTR_truncation_+_exon_loss_variant", "5_prime_UTR_variant", "initiator_codon_variant", "None", "chromosomal_deletion"
)
names(classic_snpeff_effects_to_so) <- c(
    "CDS", "CODON_CHANGE", "CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "DOWNSTREAM", "EXON", "EXON_DELETED", "FRAME_SHIFT", "GENE", "INTERGENIC", "INTERGENIC_REGION", "INTERGENIC_CONSERVED", "INTRAGENIC", "INTRON", "INTRON_CONSERVED", "NON_SYNONYMOUS_CODING", "RARE_AMINO_ACID", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "SPLICE_SITE_REGION", "START_GAINED", "START_LOST", "STOP_GAINED", "STOP_LOST", "SYNONYMOUS_CODING", "SYNONYMOUS_START", "SYNONYMOUS_STOP", "TRANSCRIPT", "UPSTREAM", "UTR_3_DELETED", "UTR_3_PRIME", "UTR_5_DELETED", "UTR_5_PRIME", "NON_SYNONYMOUS_START", "NONE", "CHROMOSOME_LARGE_DELETION"
)
# translate classic effects into SO terms leaving unknown terms (possibly SO already) as is
so_effects <- sapply(ann_final$effect, function(x) trans(x, classic_snpeff_effects_to_so))

# handle further translation of effects we care about
so_effects_translation <- list(
    "non-syn", "syn",
    "deletion", "deletion", "deletion",
    "insertion", "insertion", "frame shift",
    "stop gained", "stop lost"
)
names(so_effects_translation) <- c(
    "missense_variant", "synonymous_variant",
    "disruptive_inframe_deletion", "inframe_deletion", "chromosomal_deletion",
    "disruptive_inframe_insertion", "inframe_insertion", "frameshift_variant",
    "stop_gained", "stop_lost"
)
# translate to our simple terms turning undefined terms into '?'
simple_effects <- sapply(so_effects, function(x) trans(x, so_effects_translation, replace_missing = "?"))
# complex variant effects (those that do more than one thing) are concatenated
# with either '+' (for classic terms) or '&' (for SO terms)
simple_effects[grepl("+", so_effects, fixed = TRUE)] <- "complex"
simple_effects[grepl("&", so_effects, fixed = TRUE)] <- "complex"
simple_effects[so_effects == ""] <- "non-coding"

ann_final$effect <- simple_effects
ann_final$gene <- sub("^$", "NCR", ann_final$gene)

## automatically determine gaps for the heatmap
gap_vector <- which(!(ann_final$gene[1:length(ann_final$gene) - 1] ==  # nolint
                      ann_final$gene[2:length(ann_final$gene)]))

                                        # colormanagement
my_colors <- colorRampPalette(c("grey93", "brown", "black")) #heatmap
count <- length(unique(ann_final$gene))                     #annotations (genes)
gene_color <- c(brewer.pal(brewer_color_gene_annotation, n = count))
names(gene_color) <- unique(ann_final$gene)

                                        # colormanagement annotations (effect)
## Define the full set of colors for each effect that we can encounter
## This is not bulletproof. The effect names given here were swapped into the
## data (see above substitutions in ann_final$effect) and so are hard-coded,
## as well as their preferred colors.

all_colors <- data.frame(
    color = c("white", "green", "orange", "red",
              "black", "grey", "yellow", "blue", "purple", "brown"),
    name = c("non-coding", "syn", "non-syn", "deletion",
             "frame shift", "stop gained", "stop lost", "insertion",
             "complex", "?"))
## Reduce the full set to just those that we want
detected_effects <- unique(ann_final$effect)
subset_colors <- subset(all_colors, name %in% detected_effects)
effect_color <- subset_colors$color
names(effect_color) <- subset_colors$name
color_list <- list(gene_color = gene_color, effect_color = effect_color)
names(color_list) <- c("gene", "effect")

                                        # visualize heatmap
if (pheat_number_of_clusters > length(samples$ids)) {
    print(paste0("[INFO] Number of clusters: User-specified clusters (",
                 pheat_number_of_clusters,
                 ") is greater than the number of samples (",
                 length(samples$ids), ")"))
    pheat_number_of_clusters <- length(samples$ids)
    print(paste0("[INFO] Number of clusters: now set to ",
                 pheat_number_of_clusters))
}


                                        # Fix Labels
## Prettify names, check for label parity between final and ann_final
fix_label <- function(name) {
    ##' Reduce: 424 AGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTT A
    ##'     to: 424 AGT… > A
    cols <- unlist(str_split(name, " "))
    ## first 3 are POS REF ALT, and the rest are optional differences
    pos_ref_alt <- cols[1:3]
    rest <- ""
    if (length(cols) > 3) {
        rest <- paste0(" :: ", paste(cols[4:length(cols)], sep = " "))
    }
    ## Trim the REF or ALT if too long
    if (str_length(pos_ref_alt[2]) > 3) {
        pos_ref_alt[2] <- paste0(substring(pos_ref_alt[2], 1, 3), "…")
    }
    if (str_length(pos_ref_alt[3]) > 3) {
        pos_ref_alt[3] <- paste0(substring(pos_ref_alt[3], 1, 3), "…")
    }
    ## Join required
    new_name <- paste0(pos_ref_alt[1], " ",
                       pos_ref_alt[2], " > ",
                       pos_ref_alt[3])
    ## Join rest
    new_name <- paste0(new_name, " ", paste(rest))
}

colnames(final) <- sapply(colnames(final), fix_label)
rownames(ann_final) <- sapply(rownames(ann_final), fix_label)
## sanity test
stopifnot(all(colnames(final) %in% rownames(ann_final)))


                                        # Perform Plotting
get_plot_dims <- function(heat_map) {
    ## get the dimensions of a pheatmap object
    ## useful for plot formats that can't be written to a file directly, but
    ## for which we need to set up a plotting device
    ## source: https://stackoverflow.com/a/61876386
    plot_height <- sum(sapply(heat_map$gtable$heights,
                              grid::convertHeight, "in"))
    plot_width  <- sum(sapply(heat_map$gtable$widths,
                              grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

height <- round(max(c(max(c(
    16 * (length(unique(ann_final$effect)) +
        length(unique(ann_final$gene))), 160)) /
    nrow(final), 15)))
width <- round(ratio * height)


if (!(out_ext %in% c("svg", "jpeg", "png", "pdf"))) {
    stop("Unknown extension: ", ext, ", aborting.")
}
plot_device <- get(out_ext)


## A constant scaling factor based on the calculated dimensions
## above does not work for PNG, so we resort to feeding pheatmap
## with a direct filename
plot_filename <- NA
if (out_ext %in% c("jpeg", "png")) {
    plot_filename <- out_file
}

## SVG is not a format pheatmap knows how to write to a file directly.
## As a workaround we
## 1. create the plot object
## 2. get its dimensions
## 3. set up a svg plotting device with these dimensions
## 4. print the heatmap object to the device
hm <- pheatmap(
    final,
    color = my_colors(100),
    cellwidth = width,
    cellheight = height,
    fontsize_col = round(1 / 3 * width),
    fontsize_row = round(1 / 3 * min(c(height, width))),
    clustering_method = pheat_clustering_method,
    cluster_rows = pheat_clustering,
    cluster_cols = F,
    cutree_rows = pheat_number_of_clusters,
    annotation_col = ann_final,
    annotation_colors = color_list,
    filename = plot_filename,
    gaps_col = gap_vector
)

if (out_ext %in% c("pdf", "svg")) {
    plot_dims <- get_plot_dims(hm)
    plot_device(out_file,
                width = plot_dims$width,
                height = plot_dims$height)
    print(hm)
    dev.off()
}
