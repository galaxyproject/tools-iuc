#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(lemur)
  library(uwot)
  library(ggplot2)
})

# ---- Options ----
option_list <- list(
  make_option(c("--input_rds"), type="character", help="Input SingleCellExperiment RDS"),
  make_option(c("--design"), type="character", help="Design formula"),
  make_option(c("--contrast"), type="character", help="Contrast expression"),
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
  make_option(c("--chrom2_thresh"), type="double", default=2.5)
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Load ----
cat("Reading input:", opt$input_rds, "\n")
sce <- readRDS(opt$input_rds)

design_formula <- as.formula(opt$design)
fit <- lemur(sce, design = design_formula,
             n_embedding = opt$n_embedding,
             test_fraction = opt$test_fraction)
fit <- align_harmony(fit)
fit <- test_de(fit, contrast = eval(parse(text=opt$contrast)))

# ---- UMAP ----
umap <- uwot::umap(t(fit$embedding))
umap_df <- as_tibble(fit$colData) |>
  mutate(UMAP1 = umap[,1], UMAP2 = umap[,2])

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = patient_id, shape = condition), size = 0.5) +
  coord_fixed() + theme_minimal()
ggsave(opt$output_umap, p_umap)

# ---- Volcano ----
neighborhoods <- find_de_neighborhoods(fit, group_by = vars(patient_id, condition))
p_volcano <- neighborhoods |>
  drop_na() |>
  ggplot(aes(x = lfc, y = -log10(pval))) +
  geom_point(aes(color = adj_pval < 0.1)) +
  theme_minimal()
ggsave(opt$output_volcano, p_volcano)

# ---- Save DE without list ----
neigh_out <- neighborhoods |>
  select(-neighborhood)

write.table(as.data.frame(neigh_out), opt$output_de, sep="\t", quote=F, row.names=F)

# ---- Gene-specific UMAP + hist ----
if (!is.null(opt$sel_gene)) {
  cat("Gene-specific plots for:", opt$sel_gene, "\n")

  if (!is.null(opt$output_gene_umap)) {
    df <- tibble(umap = umap) |>
      mutate(de = assay(fit, "DE")[opt$sel_gene, ])
    p_gene_umap <- ggplot(df, aes(x = umap[,1], y = umap[,2])) +
      geom_point(aes(color = de)) +
      scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
      coord_fixed() + theme_minimal()
    ggsave(opt$output_gene_umap, p_gene_umap)
  }

  if (!is.null(opt$output_gene_hist)) {
    df <- tibble(de = assay(fit, "DE")[opt$sel_gene, ])
    p_hist <- ggplot(df, aes(x = de)) + geom_histogram(bins = 100)
    ggsave(opt$output_gene_hist, p_hist)
  }
}

# ---- Tumor block ----
if (!is.null(opt$output_chr_scatter) || !is.null(opt$output_tumor_umap) || !is.null(opt$output_tumor_neigh)) {

  if ("chromosome" %in% colnames(rowData(fit))) {
    row_data <- rowData(fit)

    chr1_expr <- colMeans(logcounts(fit)[row_data$chromosome == opt$chrom1_name, ])
    chr2_expr <- colMeans(logcounts(fit)[row_data$chromosome == opt$chrom2_name, ])

    tumor_label_df <- tibble(
      cell_id = colnames(fit),
      chrom1_expr = chr1_expr,
      chrom2_expr = chr2_expr
    ) |> mutate(
      is_tumor = chr1_expr > opt$chrom1_thresh & chrom2_expr < opt$chrom2_thresh
    )

    if (!is.null(opt$output_chr_scatter)) {
      p_chr <- ggplot(tumor_label_df, aes(x = chrom2_expr, y = chrom1_expr)) +
        geom_point(aes(color = is_tumor), size = 0.5) +
        geom_hline(yintercept = opt$chrom1_thresh) +
        geom_vline(xintercept = opt$chrom2_thresh) +
        labs(x = paste0(opt$chrom2_name, " expr"), y = paste0(opt$chrom1_name, " expr")) +
        theme_minimal()
      ggsave(opt$output_chr_scatter, p_chr)
    }

    if (!is.null(opt$output_tumor_umap)) {
      p_tumor_umap <- tibble(umap = umap) |>
        mutate(is_tumor = tumor_label_df$is_tumor) |>
        ggplot(aes(x = umap[,1], y = umap[,2])) +
        geom_point(aes(color = is_tumor), size = 0.5) +
        facet_wrap(vars(is_tumor)) +
        coord_fixed() + theme_minimal()
      ggsave(opt$output_tumor_umap, p_tumor_umap)
    }

    if (!is.null(opt$output_tumor_neigh)) {
      tumor_fit <- fit[, tumor_label_df$is_tumor]
      tum_nei <- find_de_neighborhoods(
        tumor_fit,
        group_by = vars(patient_id, condition),
        verbose = FALSE
      )
      write.table(as.data.frame(tum_nei), opt$output_tumor_neigh, sep="\t", quote=F, row.names=F)
    }

  } else {
    cat("rowData has no 'chromosome' column. Tumor plots skipped.\n")
  }
}

cat("LEMUR pipeline finished successfully.\n")
