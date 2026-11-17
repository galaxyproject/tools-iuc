#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(rlang)
    library(SingleCellExperiment)
    library(lemur)
    library(uwot)
    library(ggplot2)
})

set.seed(42)

save_plot <- function(filename, plot, format = "pdf", width = 15, height = 12, dpi = 300) {
    ggsave(filename, plot, device = tolower(format), width = width, height = height, dpi = dpi, units = "cm", bg = "white")
}

# ---- Command-line options ----
option_list <- list(
    make_option(c("--input_rds"), type = "character"),
    make_option(c("--meta_table"), type = "character"),
    make_option(c("--cell_id_column"), type = "character"),
    make_option(c("--condition_column"), type = "character"),
    make_option(c("--covariate_columns"), type = "character"),
    make_option(c("--contrast_condition"), type = "character"),
    make_option(c("--reference_condition"), type = "character"),
    make_option(c("--n_embedding"), type = "integer", default = 15),
    make_option(c("--test_fraction"), type = "double", default = 0.2),
    make_option(c("--output_umap"), type = "character"),
    make_option(c("--output_volcano"), type = "character"),
    make_option(c("--output_de"), type = "character"),
    make_option(c("--sel_gene"), type = "character", default = NULL),
    make_option(c("--output_gene_umap"), type = "character", default = NULL),
    make_option(c("--output_gene_hist"), type = "character", default = NULL),
    make_option(c("--output_chr_scatter"), type = "character", default = NULL),
    make_option(c("--output_tumor_umap"), type = "character", default = NULL),
    make_option(c("--output_tumor_neigh"), type = "character", default = NULL),
    make_option(c("--chrom1_name"), type = "character", default = "7"),
    make_option(c("--chrom2_name"), type = "character", default = "10"),
    make_option(c("--chrom1_thresh"), type = "double", default = 0.8),
    make_option(c("--chrom2_thresh"), type = "double", default = 2.5),
    make_option(c("--use_harmony"), type = "character", default = "yes"),
    make_option(c("--plot_format"), type = "character", default = "pdf"),
    make_option(c("--plot_width"), type = "double", default = 6),
    make_option(c("--plot_height"), type = "double", default = 5),
    make_option(c("--run_tumor_analysis"), type = "character", default = "no"),
    make_option(c("--tumor_annotation_column"), type = "character", default = "chromosome"),
    make_option(c("--linear_coefficient_estimator"), type = "character", default = "linear"),
    make_option(c("--use_assay"), type = "character", default = "logcounts"),
    make_option(c("--consider"), type = "character", default = "embedding+linear"),
    make_option(c("--group_by_columns"), type = "character"),
    make_option(c("--test_method"), type = "character", default = "glmGamPoi"),
    make_option(c("--n_threads"), type = "integer", default = 1)
)

opt <- parse_args(OptionParser(option_list = option_list))

sce <- readRDS(opt$input_rds)

# ---- Load and prepare metadata ----
meta <- read.delim(opt$meta_table, sep = "\t", check.names = FALSE)

opt$cell_id_column <- as.integer(opt$cell_id_column)
opt$condition_column <- as.integer(opt$condition_column)
opt$covariate_columns <- if (!is.null(opt$covariate_columns) && nzchar(opt$covariate_columns)) {
    raw_cov <- strsplit(opt$covariate_columns, ",")[[1]]
    idx <- suppressWarnings(as.integer(raw_cov))
    if (any(is.na(idx)) || any(idx < 1) || any(idx > ncol(meta))) {
        stop(sprintf(
            "Invalid covariate column index/indices: '%s'. Expected integer values between 1 and %d (the number of columns in the metadata table).",
            paste(raw_cov, collapse = ", "),
            ncol(meta)
        ))
    }
    idx
} else {
    NULL
}

# group_by is required: find_de_neighborhoods() groups cells into pseudobulk
# samples by these column(s) (typically the replication unit, e.g. patient ID).
opt$group_by_columns <- if (!is.null(opt$group_by_columns) && nzchar(opt$group_by_columns)) {
    raw_gb <- strsplit(opt$group_by_columns, ",")[[1]]
    gidx <- suppressWarnings(as.integer(raw_gb))
    if (any(is.na(gidx)) || any(gidx < 1) || any(gidx > ncol(meta))) {
        stop(sprintf(
            "Invalid group-by column index/indices: '%s'. Expected integer values between 1 and %d (the number of columns in the metadata table).",
            paste(raw_gb, collapse = ", "),
            ncol(meta)
        ))
    }
    gidx
} else {
    stop(paste0(
        "At least one group-by column is required. find_de_neighborhoods() groups ",
        "cells into pseudobulk samples by this variable (typically the replication ",
        "unit, e.g. patient or sample ID); the condition and any design covariates ",
        "are added automatically."
    ))
}

cell_id_colname <- colnames(meta)[opt$cell_id_column]
condition_name <- colnames(meta)[opt$condition_column]
batch_names <- if (!is.null(opt$covariate_columns)) colnames(meta)[opt$covariate_columns] else NULL
group_by_names <- colnames(meta)[opt$group_by_columns]

rownames(meta) <- meta[[cell_id_colname]]

if (!all(meta[[cell_id_colname]] %in% colnames(sce))) {
    stop(sprintf("Cell IDs in metadata column '%s' do not match SCE colnames", cell_id_colname))
}
stopifnot(condition_name %in% colnames(meta))
if (!is.null(batch_names)) stopifnot(all(batch_names %in% colnames(meta)))

# The subset check above only guarantees metadata cell IDs are a subset of the
# SCE columns. Also check the reverse: every SCE column must have a matching
# metadata row, otherwise meta[colnames(sce), ] would insert silent NA rows.
if (!all(colnames(sce) %in% rownames(meta))) {
    stop(sprintf(
        "Some cells in the SCE object are absent from metadata column '%s'; every SCE column must have a matching metadata row.",
        cell_id_colname
    ))
}

meta <- meta[colnames(sce), , drop = FALSE]

# Make column names syntactically valid so names with spaces do not break the
# design formula. This runs after the index-based column selection above, so it
# does not interfere with the columns the user chose by index.
colnames(meta) <- make.names(colnames(meta))
condition_name <- make.names(condition_name)
batch_names <- if (!is.null(batch_names)) make.names(batch_names) else NULL
group_by_names <- make.names(group_by_names)

colData(sce) <- S4Vectors::DataFrame(meta)

# ---- Build design formula ----
# Covariates (batch_names) may contain multiple columns; the condition goes last.
design_formula <- as.formula(
    paste("~", paste(c(batch_names, condition_name), collapse = " + "))
)

# ---- Run LEMUR ----
fit <- lemur(
    sce,
    design = design_formula,
    n_embedding = opt$n_embedding,
    test_fraction = opt$test_fraction,
    linear_coefficient_estimator = opt$linear_coefficient_estimator,
    use_assay = opt$use_assay
)

# ---- Harmony alignment (optional) ----
if (opt$use_harmony == "yes") {
    message("Running Harmony alignment...")
    fit <- align_harmony(fit)
} else {
    message("Harmony alignment skipped.")
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

fit <- test_de(fit, contrast = !!contrast_expr, consider = opt$consider)
cat("Differential expression test completed.\n")

# ---- UMAP plot ----
umap <- uwot::umap(reducedDim(fit, "embedding"), n_threads = opt$n_threads)
umap_df <- as_tibble(fit$colData) |> mutate(UMAP1 = umap[, 1], UMAP2 = umap[, 2])
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = if (!is.null(batch_names)) .data[[batch_names[1]]] else NULL, shape = .data[[condition_name]]), size = 0.5, na.rm = TRUE) +
    coord_fixed() +
    theme_minimal()
save_plot(opt$output_umap, p_umap, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)

# ---- Volcano plot ----
# group_by is only the replication unit; find_de_neighborhoods() adds the design
# formula variables (condition + covariates) to the pseudobulk grouping itself.
group_vars <- vars(!!!rlang::syms(group_by_names))
neighborhoods <- find_de_neighborhoods(fit, group_by = group_vars, test_method = opt$test_method)
if (all(c("lfc", "pval", "adj_pval") %in% colnames(neighborhoods))) {
    p_volcano <- neighborhoods |>
        drop_na() |>
        ggplot(aes(x = lfc, y = -log10(pval))) +
        geom_point(aes(color = adj_pval < 0.1)) +
        theme_minimal()
    save_plot(opt$output_volcano, p_volcano, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)
} else {
    message(
        "Skipping volcano plot: test_method = 'none' means find_de_neighborhoods() ",
        "did not compute lfc/pval/adj_pval statistics, so there is nothing to plot."
    )
}
neigh_out <- neighborhoods |> select(-neighborhood)
write.table(as.data.frame(neigh_out), opt$output_de, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Gene-specific plots ----
if (!is.null(opt$sel_gene)) {
    if (!opt$sel_gene %in% rownames(fit)) {
        stop(sprintf(
            "Selected gene '%s' was not found in the rownames of the fitted model.",
            opt$sel_gene
        ))
    }
    df <- tibble(umap = umap) |> mutate(de = assay(fit, "DE")[opt$sel_gene, ])
    if (!is.null(opt$output_gene_umap)) {
        p_gene_umap <- ggplot(df, aes(x = umap[, 1], y = umap[, 2])) +
            geom_point(aes(color = de)) +
            scale_color_gradient2(low = "blue", high = "red") +
            coord_fixed() +
            theme_minimal()
        save_plot(opt$output_gene_umap, p_gene_umap, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)
    }
    if (!is.null(opt$output_gene_hist)) {
        p_hist <- ggplot(df, aes(x = de)) +
            geom_histogram(bins = 100)
        save_plot(opt$output_gene_hist, p_hist, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)
    }
}

# ---- Tumor-specific plots (optional, GBM chromosome gain/loss example) ----
if (opt$run_tumor_analysis == "yes") {
    tumor_col <- opt$tumor_annotation_column
    if (tumor_col %in% colnames(rowData(fit))) {
        row_data <- rowData(fit)
        available_chrom_values <- unique(row_data[[tumor_col]])
        if (!(opt$chrom1_name %in% available_chrom_values)) {
            stop(sprintf(
                "Chromosome value '%s' (gained chromosome) was not found in rowData column '%s'. Available values: %s.",
                opt$chrom1_name, tumor_col, paste(available_chrom_values, collapse = ", ")
            ))
        }
        if (!(opt$chrom2_name %in% available_chrom_values)) {
            stop(sprintf(
                "Chromosome value '%s' (lost chromosome) was not found in rowData column '%s'. Available values: %s.",
                opt$chrom2_name, tumor_col, paste(available_chrom_values, collapse = ", ")
            ))
        }
        chr1_expr <- colMeans(assay(fit, opt$use_assay)[row_data[[tumor_col]] == opt$chrom1_name, ])
        chr2_expr <- colMeans(assay(fit, opt$use_assay)[row_data[[tumor_col]] == opt$chrom2_name, ])
        tumor_label_df <- tibble(cell_id = colnames(fit), chrom1_expr = chr1_expr, chrom2_expr = chr2_expr) |>
            mutate(is_tumor = chrom1_expr > opt$chrom1_thresh & chrom2_expr < opt$chrom2_thresh)

        if (!is.null(opt$output_chr_scatter)) {
            p_chr <- ggplot(tumor_label_df, aes(x = chrom2_expr, y = chrom1_expr)) +
                geom_point(aes(color = is_tumor), size = 0.5) +
                geom_hline(yintercept = opt$chrom1_thresh) +
                geom_vline(xintercept = opt$chrom2_thresh) +
                labs(x = paste0(opt$chrom2_name, " expr"), y = paste0(opt$chrom1_name, " expr")) +
                theme_minimal()
            save_plot(opt$output_chr_scatter, p_chr, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)
        }

        if (!is.null(opt$output_tumor_umap)) {
            p_tumor_umap <- tibble(umap = umap) |>
                mutate(is_tumor = tumor_label_df$is_tumor) |>
                ggplot(aes(x = umap[, 1], y = umap[, 2])) +
                geom_point(aes(color = is_tumor), size = 0.5) +
                facet_wrap(vars(is_tumor)) +
                coord_fixed() +
                theme_minimal()
            save_plot(opt$output_tumor_umap, p_tumor_umap, format = opt$plot_format, width = opt$plot_width, height = opt$plot_height)
        }

        if (!is.null(opt$output_tumor_neigh)) {
            tumor_fit <- fit[, tumor_label_df$is_tumor]
            tumor_neigh <- find_de_neighborhoods(tumor_fit, group_by = group_vars, test_method = opt$test_method)
            tumor_neigh_out <- tumor_neigh |> select(-neighborhood)
            write.table(as.data.frame(tumor_neigh_out), opt$output_tumor_neigh, sep = "\t", quote = FALSE, row.names = FALSE)
        }
    } else {
        stop(paste0(
            "The tumor annotation column '", tumor_col, "' does not exist ",
            "in the input data. Available columns in rowData: ",
            paste(colnames(rowData(fit)), collapse = ", "),
            ". Check that your input file has the correct format and contains ",
            "this column, or adjust the 'Tumor annotation column' parameter to ",
            "match your data."
        ))
    }
}

cat("LEMUR pipeline completed successfully.\n")
