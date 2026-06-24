#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(lemur)
  library(uwot)
  library(ggplot2)
})

save_plot <- function(filename, plot, format = "pdf", width = NULL, height = NULL, dpi = 300) {
  ext <- tolower(format)
  out <- sub("\\.pdf$", paste0(".", ext), filename)
  if (is.null(width)) width <- 6
  if (is.null(height)) height <- 5

  ggsave(out, plot, device = ext, width = width, height = height, dpi = dpi, units = "in", bg = "white")
}

# ---- Command-line options ----
option_list <- list(
  make_option(c("--input_rds"), type="character"),
  make_option(c("--meta_table"), type="character"),
  make_option(c("--cell_id_column"), type="character"),
  make_option(c("--condition_column"), type="character"),
  make_option(c("--batch_column"), type="character"),
  make_option(c("--contrast_condition"), type="character"),
  make_option(c("--reference_condition"), type="character"),
  make_option(c("--n_embedding"), type="integer", default=15),
  make_option(c("--test_fraction"), type="double", default=0.5),
  make_option(c("--output_umap"), type="character"),
  make_option(c("--output_volcano"), type="character"),
  make_option(c("--output_de"), type="character"),
  make_option(c("--sel_gene"), type="character", default=NULL),
  make_option(c("--output_gene_umap"), type="character", default=NULL),
  make_option(c("--output_gene_hist"), type="character", default=NULL),
  make_option(c("--output_chr_scatter"), type="character", default=NULL),
  make_option(c("--output_tumor_umap"), type="character", default=NULL),
  make_option(c("--output_tumor_neigh"), type="character", default=NULL),
  make_option(c("--chrom1_name"), type="character", default="7"),
  make_option(c("--chrom2_name"), type="character", default="10"),
  make_option(c("--chrom1_thresh"), type="double", default=0.8),
  make_option(c("--chrom2_thresh"), type="double", default=2.5),
  make_option(c("--use_harmony"), type="character", default="yes"),
  make_option(c("--plot_format"), type="character", default="pdf"),
  make_option(c("--plot_width"), type="double", default=6),
  make_option(c("--plot_height"), type="double", default=5),
  make_option(c("--tumor_annotation_column"), type="character", default="chromosome")
)

opt <- parse_args(OptionParser(option_list=option_list))

sce <- readRDS(opt$input_rds)

# ---- Load and prepare metadata ----
meta <- read.delim(opt$meta_table, sep = "\t", check.names = FALSE)

opt$cell_id_column <- as.integer(opt$cell_id_column)
opt$condition_column <- as.integer(opt$condition_column)
opt$batch_column <- if (nzchar(opt$batch_column)) as.integer(opt$batch_column) else NULL

cell_id_colname <- colnames(meta)[opt$cell_id_column]
condition_name <- colnames(meta)[opt$condition_column]
batch_name <- if (!is.null(opt$batch_column)) colnames(meta)[opt$batch_column] else NULL

rownames(meta) <- meta[[cell_id_colname]]

stopifnot(cell_id_colname %in% colnames(meta))
stopifnot(condition_name %in% colnames(meta))
if (!is.null(batch_name)) stopifnot(batch_name %in% colnames(meta))

stopifnot(all(colnames(sce) %in% rownames(meta)))
meta <- meta[!is.na(meta[[condition_name]]), , drop = FALSE]
meta <- meta[colnames(sce), , drop = FALSE]
colData(sce) <- S4Vectors::DataFrame(meta)

# ---- Build design formula ----
design_formula <- if (!is.null(batch_name)) {
  as.formula(sprintf("~ %s + %s", batch_name, condition_name))
} else {
  as.formula(sprintf("~ %s", condition_name))
}

# ---- Run LEMUR ----
fit <- lemur(
  sce,
  design = design_formula,
  n_embedding = opt$n_embedding,
  test_fraction = opt$test_fraction
)

# ---- Harmony alignment (optional) ----
if (opt$use_harmony == "yes") {
  message("Running Harmony alignment...")
  fit <- align_harmony(fit)
} else {
  message("Harmony alignment skipped.")
}

# ---- Check and fix embedding dimensions ----
embedding <- reducedDim(fit, "embedding")

if (nrow(embedding) < ncol(embedding)) {
  message("Embedding appears to have shape [cells, components]. Transposing...")
  reducedDim(fit, "embedding") <- t(embedding)
} else {
  message("Embedding shape is correct [components, cells].")
}

# ---- Differential expression contrast ----
contrast_expr <- substitute(
  cond(condition = A) - cond(condition = B),
  list(
    condition = as.name(condition_name),
    A = opt$contrast_condition,
    B = opt$reference_condition
  )
)

print(contrast_expr)
print(fit)

cat("Condition column:", condition_name, "\n")
cat("Unique values in the condition column:\n")
print(unique(colData(fit)[[condition_name]]))
cat("First values:\n")
print(head(colData(fit)[[condition_name]]))

fit <- test_de(fit, contrast = contrast_expr)
cat("Differential expression test completed.\n")

# ---- UMAP plot ----
umap <- uwot::umap(t(reducedDim(fit, "embedding")))
umap_df <- as_tibble(fit$colData) |> mutate(UMAP1 = umap[,1], UMAP2 = umap[,2])
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = if (!is.null(batch_name)) .data[[batch_name]] else NULL, shape = .data[[condition_name]]), size = 0.5, na.rm=TRUE) +
  coord_fixed() + theme_minimal()
save_plot(opt$output_umap, p_umap, format = opt$plot_format)

# ---- Volcano plot ----
group_vars <- if (!is.null(batch_name)) vars(!!sym(batch_name), !!sym(condition_name)) else vars(!!sym(condition_name))
neighborhoods <- find_de_neighborhoods(fit, group_by = group_vars)
p_volcano <- neighborhoods |> drop_na() |> ggplot(aes(x = lfc, y = -log10(pval))) +
  geom_point(aes(color = adj_pval < 0.1)) + theme_minimal()
save_plot(opt$output_volcano, p_volcano, format = opt$plot_format)
neigh_out <- neighborhoods |> select(-neighborhood)
write.table(as.data.frame(neigh_out), opt$output_de, sep="\t", quote=FALSE, row.names=FALSE)

# ---- Gene-specific plots ----
if (!is.null(opt$sel_gene)) {
  df <- tibble(umap = umap) |> mutate(de = assay(fit, "DE")[opt$sel_gene, ])
  if (!is.null(opt$output_gene_umap)) {
    p_gene_umap <- ggplot(df, aes(x = umap[,1], y = umap[,2])) +
      geom_point(aes(color = de)) +
      scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
      coord_fixed() + theme_minimal()
    save_plot(opt$output_gene_umap, p_gene_umap, format = opt$plot_format)
  }
  if (!is.null(opt$output_gene_hist)) {
    p_hist <- ggplot(df, aes(x = de)) + geom_histogram(bins = 100)
    save_plot(opt$output_gene_hist, p_hist, format = opt$plot_format)
  }
}

# ---- Tumor-specific plots ----
if (!is.null(opt$output_chr_scatter) || !is.null(opt$output_tumor_umap) || !is.null(opt$output_tumor_neigh)) {
  tumor_col <- opt$tumor_annotation_column
  if (tumor_col %in% colnames(rowData(fit))) {
    row_data <- rowData(fit)
    chr1_expr <- colMeans(logcounts(fit)[row_data[[tumor_col]] == opt$chrom1_name, ])
    chr2_expr <- colMeans(logcounts(fit)[row_data[[tumor_col]] == opt$chrom2_name, ])
    tumor_label_df <- tibble(cell_id = colnames(fit), chrom1_expr, chr2_expr) |>
      mutate(is_tumor = chrom1_expr > opt$chrom1_thresh & chr2_expr < opt$chrom2_thresh)

    if (!is.null(opt$output_chr_scatter)) {
      p_chr <- ggplot(tumor_label_df, aes(x = chrom2_expr, y = chrom1_expr)) +
        geom_point(aes(color = is_tumor), size = 0.5) +
        geom_hline(yintercept = opt$chrom1_thresh) +
        geom_vline(xintercept = opt$chrom2_thresh) +
        labs(x = paste0(opt$chrom2_name, " expr"), y = paste0(opt$chrom1_name, " expr")) +
        theme_minimal()
      save_plot(opt$output_chr_scatter, p_chr, format = opt$plot_format)
    }

    if (!is.null(opt$output_tumor_umap)) {
      p_tumor_umap <- tibble(umap = umap) |> 
        mutate(is_tumor = tumor_label_df$is_tumor) |>
        ggplot(aes(x = umap[,1], y = umap[,2])) +
        geom_point(aes(color = is_tumor), size = 0.5) +
        facet_wrap(vars(is_tumor)) +
        coord_fixed() + theme_minimal()
      save_plot(opt$output_tumor_umap, p_tumor_umap, format = opt$plot_format)
    }

    if (!is.null(opt$output_tumor_neigh)) {
      tumor_fit <- fit[, tumor_label_df$is_tumor]
      tumor_neigh <- find_de_neighborhoods(tumor_fit, group_by = group_vars)
      write.table(as.data.frame(tumor_neigh), opt$output_tumor_neigh, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  } else {
    message("The specified tumor annotation column is not present in rowData. Skipping tumor plots.")
  }
}

cat("LEMUR pipeline completed successfully.\n")
