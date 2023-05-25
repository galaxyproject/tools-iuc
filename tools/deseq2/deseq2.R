#!/usr/bin/env Rscript

# A command-line interface to DESeq2 for use with Galaxy
# written by Bjoern Gruening and modified by Michael Love 2016.03.30
#
# This argument is required:
#
#   'factors' a JSON list object from Galaxy
#
# the output file has columns:
#
#   baseMean (mean normalized count)
#   log2FoldChange (by default a moderated LFC estimate)
#   lfcSE (the standard error)
#   stat (the Wald statistic)
#   pvalue (p-value from comparison of Wald statistic to a standard Normal)
#   padj (adjusted p-value, Benjamini Hochberg correction on genes which pass the mean count filter)
#
# the first variable in 'factors' will be the primary factor.
# the levels of the primary factor are used in the order of appearance in factors.
#
# by default, levels in the order A,B,C produces a single comparison of B vs A, to a single file 'outfile'
#
# for the 'many_contrasts' flag, levels in the order A,B,C produces comparisons C vs A, B vs A, C vs B,
# to a number of files using the 'outfile' prefix: 'outfile.condition_C_vs_A' etc.
# all plots will still be sent to a single PDF, named by the arg 'plots', with extra pages.
#
# fit_type is an integer valued argument, with the options from ?estimateDisperions
#   1 "parametric"
#   2 "local"
#   3 "mean"

# setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, FALSE)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
# The line below breaks test #11 and as if not needed anymore?
# loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")


library("getopt")
library("tools")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "cores", "s", 0, "integer",
  "batch_factors", "w", 1, "character",
  "outfile", "o", 1, "character",
  "countsfile", "n", 1, "character",
  "sizefactorsfile", "F", 1, "character",
  "rlogfile", "r", 1, "character",
  "vstfile", "v", 1, "character",
  "header", "H", 0, "logical",
  "factors", "f", 1, "character",
  "files_to_labels", "l", 1, "character",
  "plots", "p", 1, "character",
  "tximport", "i", 0, "logical",
  "txtype", "y", 1, "character",
  "tx2gene", "x", 1, "character", # a space-sep tx-to-gene map or GTF/GFF3 file
  "esf", "e", 1, "character",
  "fit_type", "t", 1, "integer",
  "many_contrasts", "m", 0, "logical",
  "outlier_replace_off", "a", 0, "logical",
  "outlier_filter_off", "b", 0, "logical",
  "auto_mean_filter_off", "c", 0, "logical",
  "use_beta_priors", "d", 0, "logical",
  "alpha_ma", "A", 1, "numeric",
  "prefilter", "P", 0, "logical",
  "prefilter_value", "V", 1, "numeric"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

# enforce the following required arguments
if (is.null(opt$outfile)) {
  cat("'outfile' is required\n")
  q(status = 1)
}
if (is.null(opt$factors)) {
  cat("'factors' is required\n")
  q(status = 1)
}

verbose <- is.null(opt$quiet)

source_local <- function(fname) {
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}

source_local("get_deseq_dataset.R")

suppressPackageStartupMessages({
  library("DESeq2")
  library("RColorBrewer")
  library("gplots")
})

if (opt$cores > 1) {
  library("BiocParallel")
  register(MulticoreParam(opt$cores))
  parallel <- TRUE
} else {
  parallel <- FALSE
}

# build or read sample table

trim <- function(x) gsub("^\\s+|\\s+$", "", x)

# switch on if 'factors' was provided:
library("rjson")
parser <- newJSONParser()
parser$addData(opt$factors)
factor_list <- parser$getObject()
filenames_to_labels <- fromJSON(opt$files_to_labels)
factors <- sapply(factor_list, function(x) x[[1]])
primary_factor <- factors[1]
filenames_in <- unname(unlist(factor_list[[1]][[2]]))
labs <- unname(unlist(filenames_to_labels[basename(filenames_in)]))
sample_table <- data.frame(
  sample = basename(filenames_in),
  filename = filenames_in,
  row.names = filenames_in,
  stringsAsFactors = FALSE
)
for (factor in factor_list) {
  factor_name <- trim(factor[[1]])
  sample_table[[factor_name]] <- character(nrow(sample_table))
  lvls <- sapply(factor[[2]], function(x) names(x))
  for (i in seq_along(factor[[2]])) {
    files <- factor[[2]][[i]][[1]]
    sample_table[files, factor_name] <- trim(lvls[i])
  }
  sample_table[[factor_name]] <- factor(sample_table[[factor_name]], levels = lvls)
}
rownames(sample_table) <- labs

design_formula <- as.formula(paste("~", paste(rev(factors), collapse = " + ")))

# these are plots which are made once for each analysis
generate_generic_plots <- function(dds, factors) {
  library("ggplot2")
  library("ggrepel")
  library("pheatmap")

  rld <- rlog(dds)
  p <- plotPCA(rld, intgroup = rev(factors))
  print(p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = factor(colnames(dds))), size = 3)  + geom_point())
  dat <- assay(rld)
  dists_rl <- dist(t(dat))
  mat <- as.matrix(dists_rl)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(
    mat,
    clustering_distance_rows = dists_rl,
    clustering_distance_cols = dists_rl,
    col = colors,
    main = "Sample-to-sample distances"
  )
  plotDispEsts(dds, main = "Dispersion estimates")
}

# these are plots which can be made for each comparison, e.g.
# once for C vs A and once for B vs A
generate_specific_plots <- function(res, threshold, title_suffix) {
  use <- res$baseMean > threshold
  if (sum(!use) == 0) {
    h <- hist(res$pvalue, breaks = 0:50 / 50, plot = FALSE)
    barplot(
      height = h$counts,
      col = "powderblue",
      space = 0,
      xlab = "p-values",
      ylab = "frequency",
      main = paste("Histogram of p-values for", title_suffix)
    )
    text(x = c(0, length(h$counts)), y = 0, label = paste(c(0, 1)), adj = c(0.5, 1.7), xpd = NA)
  } else {
    h1 <- hist(res$pvalue[!use], breaks = 0:50 / 50, plot = FALSE)
    h2 <- hist(res$pvalue[use], breaks = 0:50 / 50, plot = FALSE)
    colori <- c("filtered (low count)" = "khaki", "not filtered" = "powderblue")
    barplot(
      height = rbind(h1$counts, h2$counts),
      beside = FALSE,
      col = colori,
      space = 0,
      xlab = "p-values",
      ylab = "frequency",
      main = paste("Histogram of p-values for", title_suffix)
    )
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0, 1)), adj = c(0.5, 1.7), xpd = NA)
    legend("topright", fill = rev(colori), legend = rev(names(colori)), bg = "white")
  }
    plotMA(res, main = paste("MA-plot for", title_suffix), ylim = range(res$log2FoldChange, na.rm = TRUE), alpha = opt$alpha_ma)
}

if (verbose) {
  cat(paste("primary factor:", primary_factor, "\n"))
  if (length(factors) > 1) {
    cat(paste("other factors in design:", paste(factors[-length(factors)], collapse = ","), "\n"))
  }
  cat("\n---------------------\n")
}

dds <- get_deseq_dataset(sample_table, header = opt$header, design_formula = design_formula, tximport = opt$tximport, txtype = opt$txtype, tx2gene = opt$tx2gene)
# estimate size factors for the chosen method
if (!is.null(opt$esf)) {
    dds <- estimateSizeFactors(dds, type = opt$esf)
}

# estimate size factors for each sample
# - https://support.bioconductor.org/p/97676/
if (!is.null(opt$sizefactorsfile)) {
    nm <- assays(dds)[["avgTxLength"]]
    if (!is.null(nm)) {
        ## Recommended: takes into account tximport data
        cat("\nsize factors for samples: taking tximport data into account\n")
        size_factors <- estimateSizeFactorsForMatrix(counts(dds) / nm)
    } else {
        norm_factors <- normalizationFactors(dds)
        if (!is.null(norm_factors)) {
            ## In practice, gives same results as above.
            cat("\nsize factors for samples: no tximport data, using derived normalization factors\n")
            size_factors <- estimateSizeFactorsForMatrix(norm_factors)
        } else {
            ## If we have no other information, estimate from raw.
            cat("\nsize factors for samples: no tximport data, no normalization factors, estimating from raw data\n")
            size_factors <- estimateSizeFactorsForMatrix(counts(dds))
        }
    }
    write.table(size_factors, file = opt$sizefactorsfile, sep = "\t", col.names = FALSE, quote = FALSE)
}

apply_batch_factors <- function(dds, batch_factors) {
  rownames(batch_factors) <- batch_factors$identifier
  batch_factors <- subset(batch_factors, select = -c(identifier, condition))
  dds_samples <- colnames(dds)
  batch_samples <- rownames(batch_factors)
  if (!setequal(batch_samples, dds_samples)) {
    stop("Batch factor names don't correspond to input sample names, check input files")
  }
  dds_data <- colData(dds)
  # Merge dds_data with batch_factors using indexes, which are sample names
  # Set sort to False, which maintains the order in dds_data
  reordered_batch <- merge(dds_data, batch_factors, by.x = 0, by.y = 0, sort = FALSE)
  batch_factors <- reordered_batch[, ncol(dds_data):ncol(reordered_batch)]
  for (factor in colnames(batch_factors)) {
    dds[[factor]] <- batch_factors[[factor]]
  }
  colnames(dds) <- reordered_batch[, 1]
  return(dds)
}

if (!is.null(opt$batch_factors)) {
  batch_factors <- read.table(opt$batch_factors, sep = "\t", header = TRUE)
  dds <- apply_batch_factors(dds = dds, batch_factors = batch_factors)
  batch_design <- colnames(batch_factors)[-c(1, 2)]
  design_formula <- as.formula(paste("~", paste(c(batch_design, rev(factors)), collapse = " + ")))
  design(dds) <- design_formula
}

if (verbose) {
  cat("DESeq2 run information\n\n")
  cat("sample table:\n")
  print(sample_table[, -c(1:2), drop = FALSE])
  cat("\ndesign formula:\n")
  print(design_formula)
  cat("\n\n")
  cat(paste(ncol(dds), "samples with counts over", nrow(dds), "genes\n"))
}

# minimal pre-filtering
if (!is.null(opt$prefilter)) {
    keep <- rowSums(counts(dds)) >= opt$prefilter_value
    dds <- dds[keep, ]
}

# optional outlier behavior
if (is.null(opt$outlier_replace_off)) {
  min_rep <- 7
} else {
  min_rep <- Inf
  if (verbose) cat("outlier replacement off\n")
}
if (is.null(opt$outlier_filter_off)) {
  cooks_cutoff <- TRUE
} else {
  cooks_cutoff <- FALSE
  if (verbose) cat("outlier filtering off\n")
}

# optional automatic mean filtering
if (is.null(opt$auto_mean_filter_off)) {
  independent_filtering <- TRUE
} else {
  independent_filtering <- FALSE
  if (verbose) cat("automatic filtering on the mean off\n")
}

# shrinkage of LFCs
if (is.null(opt$use_beta_priors)) {
  beta_prior <- FALSE
  if (verbose)
    cat("Applied default - beta prior off\n")
} else {
  beta_prior <- opt$use_beta_priors
}
sprintf("use_beta_prior is set to %s", beta_prior)


# dispersion fit type
if (is.null(opt$fit_type)) {
  fit_type <- "parametric"
} else {
  fit_type <- c("parametric", "local", "mean")[opt$fit_type]
}

if (verbose) cat(paste("using disperion fit type:", fit_type, "\n"))

# run the analysis
dds <- DESeq(dds, fitType = fit_type, betaPrior = beta_prior, minReplicatesForReplace = min_rep, parallel = parallel)

# create the generic plots and leave the device open
if (!is.null(opt$plots)) {
  if (verbose) cat("creating plots\n")
  pdf(opt$plots)
  generate_generic_plots(dds, factors)
}

n <- nlevels(colData(dds)[[primary_factor]])
all_levels <- levels(colData(dds)[[primary_factor]])

if (!is.null(opt$countsfile)) {
    normalized_counts <- counts(dds, normalized = TRUE)
    write.table(normalized_counts, file = opt$countsfile, sep = "\t", col.names = NA, quote = FALSE)
}

if (!is.null(opt$rlogfile)) {
    rlog_normalized <- rlogTransformation(dds)
    rlog_normalized_mat <- assay(rlog_normalized)
    write.table(rlog_normalized_mat, file = opt$rlogfile, sep = "\t", col.names = NA, quote = FALSE)
}

if (!is.null(opt$vstfile)) {
    vst_normalized <- varianceStabilizingTransformation(dds)
    vst_normalized_mat <- assay(vst_normalized)
    write.table(vst_normalized_mat, file = opt$vstfile, sep = "\t", col.names = NA, quote = FALSE)
}


if (is.null(opt$many_contrasts)) {
  # only contrast the first and second level of the primary factor
  ref <- all_levels[1]
  lvl <- all_levels[2]
  res <- results(
    dds,
    contrast = c(primary_factor, lvl, ref),
    cooksCutoff = cooks_cutoff,
    independentFiltering = independent_filtering
  )
  if (verbose) {
    cat("summary of results\n")
    cat(paste0(primary_factor, ": ", lvl, " vs ", ref, "\n"))
    print(summary(res))
  }
  res_sorted <- res[order(res$padj), ]
  out_df <- as.data.frame(res_sorted)
  out_df$geneID <- rownames(out_df)  # nolint
  out_df <- out_df[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  filename <- opt$outfile
  write.table(out_df, file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (independent_filtering) {
    threshold <- unname(attr(res, "filterThreshold"))
  } else {
    threshold <- 0
  }
  title_suffix <- paste0(primary_factor, ": ", lvl, " vs ", ref)
  if (!is.null(opt$plots)) {
    generate_specific_plots(res, threshold, title_suffix)
  }
} else {
  # rotate through the possible contrasts of the primary factor
  # write out a sorted table of results with the contrast as a suffix
  # add contrast specific plots to the device
  for (i in seq_len(n - 1)) {
    ref <- all_levels[i]
    contrast_levels <- all_levels[(i + 1):n]
    for (lvl in contrast_levels) {
      res <- results(
        dds,
        contrast = c(primary_factor, lvl, ref),
        cooksCutoff = cooks_cutoff,
        independentFiltering = independent_filtering
      )
      res_sorted <- res[order(res$padj), ]
      out_df <- as.data.frame(res_sorted)
      out_df$geneID <- rownames(out_df)  # nolint
      out_df <- out_df[, c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
      filename <- paste0(primary_factor, "_", lvl, "_vs_", ref)
      write.table(out_df, file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      if (independent_filtering) {
        threshold <- unname(attr(res, "filterThreshold"))
      } else {
        threshold <- 0
      }
      title_suffix <- paste0(primary_factor, ": ", lvl, " vs ", ref)
      if (!is.null(opt$plots)) {
        generate_specific_plots(res, threshold, title_suffix)
      }
    }
  }
}

# close the plot device
if (!is.null(opt$plots)) {
  cat("closing plot device\n")
  dev.off()
}

cat("Session information:\n\n")

sessionInfo()
