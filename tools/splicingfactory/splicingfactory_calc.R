#!/usr/bin/env Rscript

# Minimal SplicingFactory runner - simplified for clarity
suppressPackageStartupMessages({
  library(argparse)
  library(SplicingFactory)
  library(ggplot2)
  library(SummarizedExperiment)
  library(tidyr)
})

parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE)
parser$add_argument("--genes", required = TRUE)
parser$add_argument("--coldata", default = NULL)
parser$add_argument("--method", default = "laplace")
parser$add_argument("--q", type = "double", default = 2)
parser$add_argument("--norm", default = "True")
parser$add_argument("--out_div", default = "diversity.tsv")
parser$add_argument("--out_diff", default = "difference.tsv")
parser$add_argument("--plot_dir", default = NULL)
parser$add_argument("--plot_format", default = "png")
parser$add_argument("--output_summarizedExperiment_RDS", action = "store_true")
parser$add_argument("--out_se", default = "diversity_se.rds")

args <- parser$parse_args()
norm <- as.logical(args$norm)

# Read expression matrix (TSV with header; first column = transcript IDs)
tc <- read.table(args$input, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(tc) < 2) stop("Input must be a table with at least one sample column")
transcripts <- as.character(tc[[1]])
mat <- as.matrix(tc[, -1])
mode(mat) <- "numeric"
rownames(mat) <- transcripts

# Read genes mapping: prefer two-column (transcript,gene), else single-column aligned to rows
gm <- read.table(args$genes, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
if (ncol(gm) >= 2) {
  colnames(gm)[1:2] <- c("transcript", "gene")
  genes <- gm$gene[match(rownames(mat), gm$transcript)]
} else {
  genes <- as.character(gm[[1]])
  if (length(genes) != nrow(mat)) stop("One-column gene file must have same length as transcript rows")
}

# Calculate diversity
method_args <- list(x = mat, genes = genes, method = args$method, norm = norm)
if (args$method == "tsallis") method_args$q <- args$q
se <- do.call(SplicingFactory::calculate_diversity, method_args)

# Save diversity table
div_tab <- as.data.frame(assay(se))
rd <- rowData(se)
if (!is.null(rd) && "genes" %in% colnames(rd)) div_tab$Gene <- rd$genes else div_tab$Gene <- rownames(se)
div_tab <- div_tab[, c("Gene", setdiff(colnames(div_tab), "Gene"))]
write.table(div_tab, file = args$out_div, sep = "\t", quote = FALSE, row.names = FALSE)

# Optional differential analysis when coldata provided
if (!is.null(args$coldata) && file.exists(args$coldata)) {
  cd <- read.table(args$coldata, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  rownames(cd) <- as.character(cd[[1]])
  colData(se) <- S4Vectors::DataFrame(cd[colnames(se), -1, drop = FALSE])
  # determine grouping column and control group
  cd_local <- colData(se)
  if (is.null(cd_local) || ncol(cd_local) == 0) {
    warning("No metadata columns found in coldata; skipping differential analysis")
  } else {
    candidates <- intersect(c("Condition", "condition", "Group", "group"), colnames(cd_local))
    group_col <- if (length(candidates) > 0) candidates[1] else colnames(cd_local)[1]
    group_vals <- unique(as.character(cd_local[[group_col]]))
    group_vals <- group_vals[!is.na(group_vals)]
    control_value <- NULL
    if (length(group_vals) > 0) {
      preferred <- c("normal", "control", "wt", "health")
      found_pref <- group_vals[tolower(group_vals) %in% preferred]
      if (length(found_pref) > 0) control_value <- found_pref[1] else if (length(group_vals) == 2) control_value <- sort(group_vals)[1] else control_value <- group_vals[1]
    }
    diff_res <- SplicingFactory::calculate_difference(se, samples = group_col, control = control_value)
    write.table(diff_res, file = args$out_diff, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write.table(diff_res, file = args$out_diff, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Optional RDS
if (args$output_summarizedExperiment_RDS) saveRDS(se, args$out_se)

# Plotting: vignette-style overlay of Laplace (normalized) and Gini violins
if (!is.null(args$plot_dir) && dir.exists(args$plot_dir)) {
  # compute Laplace (normalized) and Gini for plotting (do not change user's selected method)
  laplace_se <- tryCatch(
    SplicingFactory::calculate_diversity(x = mat, genes = genes, method = "laplace", norm = TRUE, verbose = FALSE),
    error = function(e) NULL
  )
  gini_se <- tryCatch(
    SplicingFactory::calculate_diversity(x = mat, genes = genes, method = "gini", verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(laplace_se) && !is.null(gini_se)) {
    # prepare laplace data
    lap_df <- as.data.frame(assay(laplace_se))
    lap_df$Gene <- if (!is.null(rowData(laplace_se)) && "genes" %in% colnames(rowData(laplace_se))) rowData(laplace_se)$genes else rownames(laplace_se)
    lap_long <- pivot_longer(lap_df, -Gene, names_to = "sample", values_to = "entropy")
    lap_long$sample_type <- ifelse(grepl("_N$", lap_long$sample), "Normal", "Tumor")
    lap_long <- tidyr::drop_na(lap_long)
    lap_long$Gene2 <- paste0(lap_long$Gene, "_", lap_long$sample_type)

    # prepare gini data
    gini_df <- as.data.frame(assay(gini_se))
    gini_df$Gene <- if (!is.null(rowData(gini_se)) && "genes" %in% colnames(rowData(gini_se))) rowData(gini_se)$genes else rownames(gini_se)
    gini_long <- pivot_longer(gini_df, -Gene, names_to = "sample", values_to = "gini")
    gini_long$sample_type <- ifelse(grepl("_N$", gini_long$sample), "Normal", "Tumor")
    gini_long <- tidyr::drop_na(gini_long)
    gini_long$Gene2 <- paste0(gini_long$Gene, "_", gini_long$sample_type)

    # mean per gene/sample_type
    lap_mean <- aggregate(lap_long$entropy, by = list(lap_long$Gene2), mean)
    colnames(lap_mean) <- c("Gene2", "mean_entropy")
    lap_mean$sample_type <- ifelse(grepl("_Normal$", lap_mean$Gene2), "Normal", "Tumor")
    lap_mean$diversity <- "Normalized Laplace entropy"

    gini_mean <- aggregate(gini_long$gini, by = list(gini_long$Gene2), mean)
    colnames(gini_mean) <- c("Gene2", "mean_gini")
    gini_mean$sample_type <- ifelse(grepl("_Normal$", gini_mean$Gene2), "Normal", "Tumor")
    gini_mean$diversity <- "Gini index"

    p <- ggplot() +
      geom_violin(data = lap_mean, aes(x = sample_type, y = mean_entropy, fill = diversity), alpha = 0.6) +
      geom_violin(data = gini_mean, aes(x = sample_type, y = mean_gini, fill = diversity), alpha = 0.6) +
      scale_fill_viridis_d(name = "Diversity") +
      coord_flip() +
      theme_minimal() +
      labs(x = "Samples", y = "Diversity")

    out_base <- file.path(args$plot_dir, ifelse(tolower(args$plot_format) == "pdf", "diversity_violin.pdf", "diversity_violin.png"))
    if (tolower(args$plot_format) == "pdf") ggsave(out_base, plot = p, device = "pdf", width = 8, height = 6) else ggsave(out_base, plot = p, device = "png", width = 8, height = 6, dpi = 100)
  }
}

q(status = 0)
