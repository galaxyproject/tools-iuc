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
    "prefilter_value", "V", 1, "numeric",
    "sample_sheet_mode", "S", 0, "logical",
    "sample_sheet", "g", 1, "character",
    "factor_columns", "j", 1, "character",
    "reference_level", "R", 1, "character",
    "target_level", "T", 1, "character",
    "collection_files", "C", 1, "character",
    "custom_design_formula", "D", 0, "logical",
    "design_formula", "G", 1, "character",
    "contrast_definition", "K", 1, "character"
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
if (is.null(opt$factors) && is.null(opt$sample_sheet_mode)) {
    cat("'factors' is required when not using sample_sheet_mode\n")
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

# Helper function to decode base64-encoded JSON
decode_base64_json <- function(encoded_str) {
    decoded_bytes <- base64enc::base64decode(encoded_str)
    decoded_str <- rawToChar(decoded_bytes)
    return(decoded_str)
}

# switch on if 'factors' was provided:
library("rjson")

if (!is.null(opt$sample_sheet_mode)) {
    # Sample sheet mode: build factor_list from sample sheet
    filenames_to_labels <- fromJSON(decode_base64_json(opt$files_to_labels))

    # Read sample sheet
    sample_sheet <- read.table(opt$sample_sheet, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    # Parse collection files
    collection_files <- strsplit(opt$collection_files, ",")[[1]]

    # First column of sample sheet should contain sample identifiers
    sample_id_col <- colnames(sample_sheet)[1]

    if (!is.null(opt$custom_design_formula)) {
        # Custom design formula mode
        # Extract variable names from formula (remove ~, +, *, :, whitespace)
        formula_str <- gsub("^~\\s*", "", opt$design_formula)
        formula_vars <- unique(trimws(unlist(strsplit(formula_str, "[+*:]"))))

        # Validate all variables exist in sample sheet
        missing_vars <- setdiff(formula_vars, colnames(sample_sheet))
        if (length(missing_vars) > 0) {
            cat(paste0("Error: Variables not found in sample sheet: ", paste(missing_vars, collapse = ", "), "\n"))
            cat(paste0("Available columns: ", paste(colnames(sample_sheet), collapse = ", "), "\n"))
            q(status = 1)
        }

        # Use all formula variables as factor columns
        factor_col_names <- formula_vars
    } else {
        # Automatic mode: Parse factor columns (comma-separated column numbers, 1-indexed)
        factor_col_nums <- as.integer(strsplit(opt$factor_columns, ",")[[1]])
        factor_col_names <- colnames(sample_sheet)[factor_col_nums]
    }

    # Validate sample sheet matches collection before building factor_list
    # Get element identifiers from collection
    collection_element_ids <- character(length(collection_files))
    for (idx in seq_along(collection_files)) {
        file <- collection_files[idx]
        file_basename <- basename(file)
        if (file_basename %in% names(filenames_to_labels)) {
            collection_element_ids[idx] <- filenames_to_labels[[file_basename]]
        } else {
            cat("Error: Sample sheet validation failed\n")
            cat(paste0("Collection file '", file_basename, "' does not have a corresponding element identifier.\n"))
            cat("This is an internal error - please report this issue.\n")
            q(status = 1)
        }
    }

    # Get sample identifiers from sample sheet
    sample_sheet_ids <- sample_sheet[[sample_id_col]]

    # Check for mismatches
    collection_not_in_sheet <- setdiff(collection_element_ids, sample_sheet_ids)
    sheet_not_in_collection <- setdiff(sample_sheet_ids, collection_element_ids)

    if (length(collection_not_in_sheet) > 0 || length(sheet_not_in_collection) > 0) {
        cat("Error: Sample sheet does not match the input collection\n\n")

        if (length(collection_not_in_sheet) > 0) {
            cat("The following samples are in the collection but NOT in the sample sheet:\n")
            for (id in collection_not_in_sheet) {
                cat(paste0("  - ", id, "\n"))
            }
            cat("\n")
        }

        if (length(sheet_not_in_collection) > 0) {
            cat("The following samples are in the sample sheet but NOT in the collection:\n")
            for (id in sheet_not_in_collection) {
                cat(paste0("  - ", id, "\n"))
            }
            cat("\n")
        }

        cat("Please ensure that:\n")
        cat(paste0("1. The first column ('", sample_id_col, "') of the sample sheet contains sample identifiers\n"))
        cat("2. These identifiers exactly match the element identifiers in your collection\n")
        cat("3. All samples in the collection are listed in the sample sheet\n")
        cat("4. All samples in the sample sheet exist in the collection\n")
        q(status = 1)
    }

    # Determine which factor will be the primary factor for contrasts
    # In custom mode, the primary factor is the last one in the formula (rightmost term)
    # In automatic mode, the first selected factor is primary (formula gets reversed)
    if (!is.null(opt$custom_design_formula)) {
        primary_factor_index <- length(factor_col_names)
    } else {
        primary_factor_index <- 1
    }

    # Build factor_list structure
    factor_list <- list()
    for (i in seq_along(factor_col_names)) {
        factor_name <- factor_col_names[i]

        # Group files by factor level, preserving the order of first appearance
        level_to_files <- list()
        level_order <- character(0) # Track order of first appearance
        for (file in collection_files) {
            element_id <- filenames_to_labels[[basename(file)]]
            # Find matching row in sample sheet
            matching_row <- which(sample_sheet[[sample_id_col]] == element_id)
            if (length(matching_row) > 0) {
                # Convert level to character to ensure consistent list indexing
                level <- as.character(sample_sheet[[factor_name]][matching_row[1]])
                if (!(level %in% names(level_to_files))) {
                    level_to_files[[level]] <- character(0)
                    level_order <- c(level_order, level) # Record order of first appearance
                }
                level_to_files[[level]] <- c(level_to_files[[level]], file)
            }
        }

        # Get all levels in order of first appearance (not alphabetically sorted)
        all_levels <- level_order

        # Handle reference and target levels for the primary contrast factor
        # In custom mode: this is the LAST factor in the formula
        # In automatic mode: this is the FIRST factor selected (formula will be reversed)
        if (i == primary_factor_index) {
            # This is the primary factor for contrasts: handle reference/target
            # Determine reference level
            if (!is.null(opt$reference_level) && nchar(opt$reference_level) > 0) {
                ref_level <- trim(opt$reference_level)
                # Validate that reference level exists
                if (!(ref_level %in% all_levels)) {
                    cat(paste0("Error: Reference level '", ref_level, "' not found in factor '", factor_name, "'. Available levels: ", paste(all_levels, collapse = ", "), "\n"))
                    q(status = 1)
                }
            } else {
                # Default: use first level encountered in sample sheet as reference
                ref_level <- all_levels[1]
            }

            # Build factor levels with reference first, then others
            # This applies to both automatic and custom modes
            if (!is.null(opt$target_level) && nchar(opt$target_level) > 0) {
                target_level <- trim(opt$target_level)
                # Validate that target level exists
                if (!(target_level %in% all_levels)) {
                    cat(paste0("Error: Target level '", target_level, "' not found in factor '", factor_name, "'. Available levels: ", paste(all_levels, collapse = ", "), "\n"))
                    q(status = 1)
                }
                # Explicit target: reference first, target second
                levels_to_use <- c(ref_level, target_level)
            } else {
                # Reference first, then others in order of appearance
                other_levels <- all_levels[all_levels != ref_level]
                levels_to_use <- c(ref_level, other_levels)
                # Target is the first non-reference level
                target_level <- other_levels[1]
            }

            # Store reference and target levels for later contrast specification
            # This ensures we use the intended levels, not whatever order files happened to be in
            primary_ref_level <- ref_level
            primary_target_level <- target_level
        } else {
            # For secondary factors, use all levels in order of appearance
            levels_to_use <- all_levels
        }

        # Build factor structure in the order specified by levels_to_use
        # Following standard DESeq2 convention: reference level is first
        factor_levels <- list()
        for (level in levels_to_use) {
            if (level %in% names(level_to_files)) {
                level_entry <- list()
                level_entry[[level]] <- level_to_files[[level]]
                factor_levels[[length(factor_levels) + 1]] <- level_entry
            }
        }

        factor_list[[length(factor_list) + 1]] <- list(factor_name, factor_levels)
    }

    # Set primary_factor for sample sheet mode
    # This is the factor that will be used for contrasts
    if (!is.null(opt$custom_design_formula)) {
        # Custom mode: primary factor is the last one in the formula (rightmost term)
        primary_factor <- factor_col_names[length(factor_col_names)]
    } else {
        # Automatic mode: primary factor is the first selected factor (formula will be reversed)
        primary_factor <- factor_col_names[1]
    }

    # Parse contrast_definition if in custom mode
    if (!is.null(opt$custom_design_formula) && !is.null(opt$contrast_definition) && nchar(opt$contrast_definition) > 0) {
        contrast_parts <- strsplit(opt$contrast_definition, ",")[[1]]
        if (length(contrast_parts) != 3) {
            cat("Error: contrast_definition must be in format 'factor,target,reference'\n")
            q(status = 1)
        }
        custom_contrast <- list(
            factor = trim(contrast_parts[1]),
            target = trim(contrast_parts[2]),
            reference = trim(contrast_parts[3])
        )
        # Validate the contrast
        if (!(custom_contrast$factor %in% factor_col_names)) {
            cat(paste0("Error: Contrast factor '", custom_contrast$factor, "' not found in design formula.\n"))
            cat(paste0("Available factors: ", paste(factor_col_names, collapse = ", "), "\n"))
            q(status = 1)
        }
    }
} else {
    # Original mode: factors provided directly
    parser <- newJSONParser()
    parser$addData(decode_base64_json(opt$factors))
    factor_list <- parser$getObject()
    filenames_to_labels <- fromJSON(decode_base64_json(opt$files_to_labels))

    # For original mode, extract reference and target levels from the first factor
    # In original mode: ref=level1 (denominator), target=level2 (numerator) -> log2(level2/level1)
    # So we swap the naming to match the contrast direction used below
    primary_factor_data <- factor_list[[1]]
    primary_levels <- sapply(primary_factor_data[[2]], function(x) names(x))
    primary_target_level <- primary_levels[1] # First level becomes "target" for the swap below
    primary_ref_level <- if (length(primary_levels) >= 2) primary_levels[2] else primary_levels[1] # Second level becomes "ref"
}

factors <- sapply(factor_list, function(x) x[[1]])
# For original mode (not sample_sheet_mode), set primary_factor to first factor
# In sample_sheet_mode, it's already set correctly above
if (is.null(opt$sample_sheet_mode)) {
    primary_factor <- factors[1]
}
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
    # Explicitly set the reference level using relevel() for the primary factor in sample_sheet_mode
    # This ensures DESeq2 knows the reference regardless of factor level order
    # In original mode, we don't relevel because the order from factor_list is already correct
    # Note: primary_factor is already set correctly in sample_sheet_mode before we get here
    if (exists("primary_factor") && factor_name == primary_factor && exists("primary_ref_level") && !is.null(opt$sample_sheet_mode)) {
        sample_table[[factor_name]] <- relevel(sample_table[[factor_name]], ref = primary_ref_level)
    }
}
rownames(sample_table) <- labs

# Build design formula
if (!is.null(opt$custom_design_formula)) {
    # Custom mode: use user-provided formula
    design_formula <- as.formula(opt$design_formula)
    # In original mode (not sample_sheet), set primary factor to last variable in formula
    # In sample_sheet_mode, it's already set correctly above
    if (is.null(opt$sample_sheet_mode)) {
        primary_factor <- factors[length(factors)]
    }
} else {
    # Automatic mode: build from selected factors
    design_formula <- as.formula(paste("~", paste(rev(factors), collapse = " + ")))
}

# these are plots which are made once for each analysis
generate_generic_plots <- function(dds, factors) {
    library("ggplot2")
    library("ggrepel")
    library("pheatmap")

    rld <- rlog(dds)
    p <- plotPCA(rld, intgroup = rev(factors))
    print(p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = factor(colnames(dds))), size = 3) + geom_point())
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

# use/estimate size factors with the chosen method
if (!is.null(opt$esf)) {
    if (opt$esf %in% list("ratio", "poscounts", "iterate")) {
        cat("Calculating size factors de novo\n")
        dds <- estimateSizeFactors(dds, type = opt$esf)
    } else {
        sf_table <- read.table(opt$esf)
        # Sort the provided size factors just in case the order differs from the input file order.
        merged_table <- merge(sample_table, sf_table, by.x = 0, by.y = 1, sort = FALSE)
        sf_values <- as.numeric(unlist(merged_table[5]))
        "sizeFactors"(dds) <- sf_values

        cat("Using user-provided size factors:\n")
        print(sf_values)
    }
} else {
    cat("No size factor was used\n")
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

    if (!is.null(opt$custom_design_formula)) {
        # Custom mode: prepend batch factors to user's formula
        user_formula_rhs <- gsub("^~\\s*", "", opt$design_formula)
        design_formula <- as.formula(paste("~", paste(c(batch_design, user_formula_rhs), collapse = " + ")))
    } else {
        # Automatic mode: prepend batch factors to generated formula
        design_formula <- as.formula(paste("~", paste(c(batch_design, rev(factors)), collapse = " + ")))
    }
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
    if (verbose) {
        cat("Applied default - beta prior off\n")
    }
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
    # Single contrast mode
    if (!is.null(opt$custom_design_formula) && exists("custom_contrast")) {
        # Use explicit custom contrast
        ref <- custom_contrast$reference
        lvl <- custom_contrast$target
        contrast_factor <- custom_contrast$factor
    } else {
        # Default: contrast using stored reference and target levels
        # This ensures we use the intended reference, not whatever order files happened to be in

        # Check if we have more than 2 levels without explicit target
        if (length(all_levels) > 2 && is.null(opt$target_level)) {
            cat("Error: Multiple factor levels detected without explicit target_level specification.\n")
            cat(paste0("Factor '", primary_factor, "' has ", length(all_levels), " levels: ", paste(all_levels, collapse = ", "), "\n"))
            cat("Please either:\n")
            cat("  1. Specify target_level to indicate which level to compare against the reference, or\n")
            cat("  2. Use many_contrasts mode to compare all levels pairwise\n")
            q(status = 1)
        }

        # Use stored reference and target levels
        # contrast = c(factor, level1, level2) computes log2(level1 / level2)
        # In original mode: variable names were swapped for backward compatibility
        # In sample_sheet_mode: use them directly (ref=denominator, target=numerator)
        if (!is.null(opt$sample_sheet_mode)) {
            # Sample sheet mode: use natural order
            ref <- primary_ref_level # Reference is denominator
            lvl <- primary_target_level # Target is numerator
        } else {
            # Original mode: swap for backward compatibility
            ref <- primary_target_level # This becomes the denominator
            lvl <- primary_ref_level # This becomes the numerator
        }
        contrast_factor <- primary_factor
    }

    res <- results(
        dds,
        contrast = c(contrast_factor, lvl, ref),
        cooksCutoff = cooks_cutoff,
        independentFiltering = independent_filtering
    )
    if (verbose) {
        cat("summary of results\n")
        cat(paste0(contrast_factor, ": ", lvl, " vs ", ref, "\n"))
        print(summary(res))
    }
    res_sorted <- res[order(res$padj), ]
    out_df <- as.data.frame(res_sorted)
    out_df$geneID <- rownames(out_df) # nolint
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
            out_df$geneID <- rownames(out_df) # nolint
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
