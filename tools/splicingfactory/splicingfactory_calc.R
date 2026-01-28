#!/usr/bin/env Rscript

library("argparse", quietly = TRUE, warn.conflicts = FALSE)
library("SplicingFactory",
  quietly = TRUE,
  warn.conflicts = FALSE
)
library("ggplot2", quietly = TRUE, warn.conflicts = FALSE)
library("SummarizedExperiment",
  quietly = TRUE,
  warn.conflicts = FALSE
)
library("tidyr", quietly = TRUE, warn.conflicts = FALSE)

options(
  show.error.messages = TRUE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)

suppressWarnings(Sys.setlocale("LC_MESSAGES", "en_US.UTF-8"))

parser <- ArgumentParser(
  description = "Galaxy runner for SplicingFactory::calculate_diversity"
)

parser$add_argument(
  "--input",
  help = "Input file: TSV matrix or RDS (tximport/SE)",
  required = TRUE
)

parser$add_argument(
  "--genes",
  help = paste0(
    "Two-column TSV mapping transcripts to genes (transcript\\tgene). ",
    "Required for TSV matrix input."
  ),
  default = NULL,
  required = FALSE
)

parser$add_argument(
  "--method",
  help = "Method: laplace, naive, gini, simpson, invsimpson, tsallis",
  default = "laplace"
)

parser$add_argument(
  "--q",
  help = "Tsallis entropy parameter q (default: 2)",
  type = "double",
  default = 2,
  required = FALSE
)

parser$add_argument(
  "--norm",
  help = "Normalize entropy: True/False",
  default = "True"
)

parser$add_argument(
  "--tpm",
  help = "Use TPM values from tximport: True/False",
  default = "False"
)

parser$add_argument(
  "--verbose",
  help = "Verbose logging",
  action = "store_true",
  default = FALSE
)

parser$add_argument(
  "--plot_dir",
  help = "Directory to save plots",
  default = NULL,
  required = FALSE
)

parser$add_argument(
  "--plot_format",
  help = "Image format for plots: png or pdf",
  default = "png",
  required = FALSE
)

parser$add_argument(
  "--diff_method",
  help = "Method to summarize per-group diversity: mean or median",
  default = "mean",
  required = FALSE
)

parser$add_argument(
  "--test",
  help = "Test for significance: wilcoxon or label_shuffling",
  default = "wilcoxon",
  required = FALSE
)

parser$add_argument(
  "--randomizations",
  help = "Number of randomizations for label shuffling test",
  type = "integer",
  default = 100,
  required = FALSE
)

parser$add_argument(
  "--pcorr",
  help = "P-value correction method (e.g. BH)",
  default = "BH",
  required = FALSE
)

parser$add_argument(
  "--out_diff",
  help = "Output differential analysis TSV file",
  default = "difference.tsv",
  required = FALSE
)


parser$add_argument(
  "--coldata",
  help = "Optional sample metadata TSV file (first column: sample IDs).",
  default = NULL,
  required = FALSE
)

parser$add_argument(
  "--out_div",
  help = "Output diversity TSV file",
  default = "diversity.tsv",
  required = FALSE
)

parser$add_argument(
  "--out_se",
  help = "Output SummarizedExperiment RDS file",
  default = "diversity_se.rds",
  required = FALSE
)

parser$add_argument(
  "--session_out",
  help = "sessionInfo() output file",
  default = "session.txt",
  required = FALSE
)

# Allow Galaxy to request whether to output the SummarizedExperiment RDS
parser$add_argument(
  "--output_summarizedExperiment_RDS",
  help = "If set, save SummarizedExperiment RDS output",
  action = "store_true",
  default = FALSE
)

cmd <- commandArgs(trailingOnly = TRUE)
raw_args <- paste(shQuote(cmd), collapse = " ")
cat("RAW_CMD_ARGS:", raw_args, "\n", file = stderr())

args <- parser$parse_args()

# Convert string booleans to R logical
norm <- as.logical(args$norm)
use_tpm <- as.logical(args$tpm)

if (args$verbose) {
  cat("Parameters:\n", file = stderr())
  cat("  Input:", args$input, "\n", file = stderr())
  cat("  Method:", args$method, "\n", file = stderr())
  cat("  Normalize:", norm, "\n", file = stderr())
  cat("  Use TPM:", use_tpm, "\n", file = stderr())
}

# ============================================================================
# Read input data
# ============================================================================

if (args$verbose) {
  cat("Reading input data...\n", file = stderr())
}


tc_df <- read.table(
  args$input,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = '"',
  check.names = FALSE,
  fill = TRUE
)
if (ncol(tc_df) < 2) {
  stop("Input matrix must contain at least one ID column and one sample column")
}
# Detect whether the first column is a transcript ID column or if the
# file is a pure numeric matrix (no ID column). If the first column is
# numeric for all rows, assume there's no ID column and use the full
# table as numeric matrix. Otherwise treat first column as transcript IDs.
first_col <- tc_df[[1]]
suppressWarnings(
  first_col_numeric <- !any(is.na(as.numeric(as.character(first_col))))
)

if (first_col_numeric) {
  # No explicit transcript ID column; use entire table as matrix
  transcript_ids <- NULL
  mat <- tc_df[, , drop = FALSE]
} else {
  transcript_ids <- as.character(first_col)
  if (any(duplicated(transcript_ids))) {
    warning("Duplicate transcript IDs found in input; keeping duplicate IDs as provided")
  }
  mat <- tc_df[, -1, drop = FALSE]
}
# Ensure numeric conversion for all sample columns
mat <- as.matrix(mat)
mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
if (!is.null(transcript_ids)) {
  rownames(mat) <- transcript_ids
}
readcounts <- mat

# Read gene mapping
# Read gene mapping: accept either two-column (transcript \t gene)
# or one-column (gene per transcript, in same order as the matrix).
gene_map_raw <- read.table(
  args$genes,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = "\""
)

if (ncol(gene_map_raw) == 1) {
  # Single-column file: treat as gene vector aligned to rows of the matrix.
  # Be tolerant of common formatting issues: drop empty/NA rows and detect
  # a possible header row if it causes an off-by-one length.
  gene_vec <- as.character(gene_map_raw[[1]])
  # Remove empty lines and NAs
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]

  # If there's exactly one extra entry, it may be a header row; drop it.
  if (length(gene_vec) == nrow(readcounts) + 1) {
    warning("One extra line in gene mapping: assuming header, dropping first row")
    gene_vec <- gene_vec[-1]
  }

  if (length(gene_vec) != nrow(readcounts)) {
    msg <- paste0(
      "Gene mapping file has 1 column but length (", length(gene_vec),
      ") does not match number of transcripts (", nrow(readcounts),
      "). Ensure the gene file has one entry per transcript row, or provide a two-column mapping."
    )
    stop(msg)
  }

  genes <- gene_vec
} else {
  # Two-or-more-column file:
  # use first two columns as transcript -> gene mapping
  gene_map <- gene_map_raw
  colnames(gene_map)[1:2] <- c("transcript", "gene")
  gene_map$transcript <- as.character(gene_map$transcript)
  gene_map$gene <- as.character(gene_map$gene)

  # Ensure order matches (transcripts must be present in mapping)
  if (!all(rownames(readcounts) %in% gene_map$transcript)) {
    stop(
      "Some transcripts in count matrix not found in",
      " gene mapping file"
    )
  }

  # Map transcripts to genes (keep transcript order)
  genes <- gene_map$gene[match(rownames(readcounts), gene_map$transcript)]
}


if (args$verbose) {
  cat(
    "  Loaded",
    nrow(readcounts),
    "transcripts,",
    ncol(readcounts),
    "samples\n",
    file = stderr()
  )
}

# ============================================================================
# Basic filtering
# ============================================================================

if (args$verbose) {
  cat("Filtering low-abundance transcripts...\n", file = stderr())
}

keep_idx <- rowSums(readcounts > 5) > 5
readcounts <- readcounts[keep_idx, ]
genes <- genes[keep_idx]

if (args$verbose) {
  cat("  Retained",
    nrow(readcounts),
    "transcripts after filtering\n",
    file = stderr()
  )
}

# ============================================================================
# Calculate diversity
# ============================================================================

if (args$verbose) {
  cat("Calculating diversity using method:",
    args$method,
    "\n",
    file = stderr()
  )
}

# Prepare method parameters
method_args <- list(
  x = readcounts,
  genes = genes,
  method = args$method,
  norm = norm,
  verbose = args$verbose
)

# Add q parameter for Tsallis
if (args$method == "tsallis") {
  method_args$q <- args$q
  if (args$verbose) {
    cat("  Using q =", args$q, "for Tsallis entropy\n", file = stderr())
  }
}

# Call calculate_diversity
result_se <- do.call(SplicingFactory::calculate_diversity, method_args)

if (args$verbose) {
  cat("Diversity calculation complete\n", file = stderr())
}

# If user provided a coldata TSV, read and attach it to colData(result_se)
if (!is.null(args$coldata) && file.exists(args$coldata)) {
  if (args$verbose) cat("Reading coldata from:", args$coldata, "\n", file = stderr())
  cd_df <- read.table(
    args$coldata,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = '"',
    check.names = FALSE
  )
  if (ncol(cd_df) < 1) {
    warning("coldata file appears empty or malformed; skipping coldata attachment")
  } else {
    sample_ids <- as.character(cd_df[[1]])
    # set rownames to sample IDs and drop the first column
    rownames(cd_df) <- sample_ids
    cd_meta <- if (ncol(cd_df) > 1) cd_df[, -1, drop = FALSE] else data.frame(row.names = sample_ids)

    # Ensure samples match columns of result_se
    se_samples <- colnames(result_se)
    if (!all(se_samples %in% rownames(cd_meta))) {
      missing_samples <- setdiff(se_samples, rownames(cd_meta))
      warning("The following samples are missing from coldata and will have NA metadata: ", paste(missing_samples, collapse = ", "))
      # create empty rows for missing samples
      for (ms in missing_samples) {
        cd_meta[ms, ] <- NA
      }
    }

    # Reorder to match result_se columns
    cd_meta_ordered <- cd_meta[se_samples, , drop = FALSE]
    # Attach as colData
    S4Vectors::metadata(result_se)$coldata_source <- args$coldata
    colData(result_se) <- S4Vectors::DataFrame(cd_meta_ordered)
    if (args$verbose) cat("Attached coldata to result_se\n", file = stderr())
  }
} else if (!is.null(args$coldata) && !file.exists(args$coldata)) {
  warning(paste0("coldata file not found: ", args$coldata))
}

# ============================================================================
# Generate plots
# ============================================================================

if (!is.null(args$plot_dir) && dir.exists(args$plot_dir)) {

  diversity_values <- assay(result_se)

  # Get gene IDs - handle different possible structures from SplicingFactory
  rd <- rowData(result_se)
  if (!is.null(rd) && "genes" %in% colnames(rd)) {
    gene_ids <- rd$genes
  } else if (!is.null(rd) && ncol(rd) > 0) {
    gene_ids <- rd[[1]]
  } else {
    gene_ids <- rownames(diversity_values)
  }

  # Create data frame for plotting
  plot_data <- as.data.frame(diversity_values)
  plot_data$Gene <- gene_ids
  plot_data_long <- pivot_longer(plot_data, -Gene, names_to = "Sample", values_to = "Diversity")

  # Auto-detect grouping column in colData (prefer 'Condition' or 'Group')
  group_col <- NULL
  cd <- colData(result_se)
  if (!is.null(cd) && ncol(cd) > 0) {
    candidates <- c("Condition", "condition", "Group", "group")
    found <- intersect(candidates, colnames(cd))
    if (length(found) > 0) {
      group_col <- found[1]
    } else {
      group_col <- colnames(cd)[1]
    }
    if (args$verbose) cat("Using grouping column:", group_col, "\n", file = stderr())

    # Attach Group vector to plotting data (align by sample names)
    if (group_col %in% colnames(cd)) {
      group_vec <- as.character(cd[[group_col]])
      plot_data_long$Group <- group_vec[match(plot_data_long$Sample, colnames(result_se))]
    }
  }

  # Build violin plot to match vignette: overlay Laplace and Gini per-sample-type
  # Compute both Laplace and Gini entropies (normalized) for plotting
  laplace_se <- SplicingFactory::calculate_diversity(readcounts, genes, method = "laplace", norm = TRUE, verbose = FALSE)
  gini_se <- SplicingFactory::calculate_diversity(readcounts, genes, method = "gini", verbose = FALSE)

  # Construct data.frames and reshape like the vignette
  laplace_data <- cbind(assay(laplace_se), Gene = rowData(laplace_se)$genes)
  laplace_data <- pivot_longer(as.data.frame(laplace_data), -Gene, names_to = "sample", values_to = "entropy")
  laplace_data$sample_type <- apply(laplace_data[, "sample", drop = FALSE], 1, function(x) ifelse(grepl("_N$", x), "Normal", "Tumor"))
  laplace_data <- tidyr::drop_na(laplace_data)
  laplace_data$Gene <- paste0(laplace_data$Gene, "_", laplace_data$sample_type)
  laplace_data$diversity <- "Normalized Laplace entropy"

  gini_data <- cbind(assay(gini_se), Gene = rowData(gini_se)$genes)
  gini_data <- pivot_longer(as.data.frame(gini_data), -Gene, names_to = "sample", values_to = "gini")
  gini_data$sample_type <- apply(gini_data[, "sample", drop = FALSE], 1, function(x) ifelse(grepl("_N$", x), "Normal", "Tumor"))
  gini_data <- tidyr::drop_na(gini_data)
  gini_data$Gene <- paste0(gini_data$Gene, "_", gini_data$sample_type)
  gini_data$diversity <- "Gini index"

  # Mean per gene/sample_type (genes encoded with sample_type suffix)
  laplace_entropy_mean <- aggregate(laplace_data$entropy, by = list(laplace_data$Gene), mean)
  colnames(laplace_entropy_mean)[2] <- "mean_entropy"
  laplace_entropy_mean <- as.data.frame(laplace_entropy_mean)
  laplace_entropy_mean$sample_type <- ifelse(grepl("_Normal$", laplace_entropy_mean[,1]), "Normal", "Tumor")
  laplace_entropy_mean$diversity <- "Normalized Laplace entropy"

  gini_mean <- aggregate(gini_data$gini, by = list(gini_data$Gene), mean)
  colnames(gini_mean)[2] <- "mean_gini"
  gini_mean <- as.data.frame(gini_mean)
  gini_mean$sample_type <- ifelse(grepl("_Normal$", gini_mean[,1]), "Normal", "Tumor")
  gini_mean$diversity <- "Gini index"

  p_violin <- ggplot() +
    geom_violin(data = laplace_entropy_mean, aes(x = sample_type, y = mean_entropy, fill = diversity), alpha = 0.6) +
    geom_violin(data = gini_mean, aes(x = sample_type, y = mean_gini, fill = diversity), alpha = 0.6) +
    scale_fill_viridis_d(name = "Diversity") +
    coord_flip() +
    theme_minimal() +
    labs(x = "Samples", y = "Diversity")

  # Save violin plot
  out_fname_base <- ifelse(args$plot_format == "pdf", "diversity_violin.pdf", "diversity_violin.png")
  out_path <- file.path(args$plot_dir, out_fname_base)
  if (args$plot_format == "pdf") {
    ggsave(out_path, plot = p_violin, device = "pdf", width = 8, height = 6)
  } else {
    ggsave(out_path, plot = p_violin, device = "png", width = 8, height = 6, dpi = 100)
  }
}

# ============================================================================
# Output diversity table
# ============================================================================

if (args$verbose) {
  cat("Writing output files...\n", file = stderr())
}

# Extract diversity values and gene info
diversity_table <- as.data.frame(assay(result_se))

# Add gene IDs - handle different possible structures
rd <- rowData(result_se)
if (!is.null(rd) && "genes" %in% colnames(rd)) {
  diversity_table$Gene <- rd$genes
} else if (!is.null(rd) && ncol(rd) > 0) {
  diversity_table$Gene <- rd[[1]]
} else {
  diversity_table$Gene <- rownames(result_se)
}

# Reorder columns to put Gene first
diversity_table <- diversity_table[
  , c("Gene", setdiff(colnames(diversity_table), "Gene"))
]

write.table(
  diversity_table,
  file = args$out_div,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (args$verbose) {
  cat("  Saved diversity table to:", args$out_div, "\n", file = stderr())
}

## Run differential analysis if a valid coldata file was provided
if (!is.null(args$coldata) && file.exists(args$coldata)) {
  if (args$verbose) cat("Running differential analysis (calculate_difference) using coldata...\n", file = stderr())

  cd <- colData(result_se)
  if (is.null(cd) || ncol(cd) == 0) {
    warning("colData not attached to SummarizedExperiment; skipping differential analysis")
  } else {
    # Auto-detect grouping column in colData (prefer 'Condition' or 'Group')
    candidates <- c("Condition", "condition", "Group", "group")
    found <- intersect(candidates, colnames(cd))
    if (length(found) > 0) {
      group_col <- found[1]
    } else {
      group_col <- colnames(cd)[1]
    }

    if (args$verbose) cat("Using samples column for calculate_difference:", group_col, "\n", file = stderr())

    # Auto-detect control group from colData values
    group_vals <- unique(as.character(cd[[group_col]]))
    group_vals <- group_vals[!is.na(group_vals)]
    control_value <- NULL
    if (length(group_vals) == 0) {
      warning("No non-NA groups found in samples column; skipping differential analysis")
    } else {
      # prefer common control labels
      preferred <- c("normal", "control", "wt", "health")
      found_pref <- group_vals[tolower(group_vals) %in% preferred]
      if (length(found_pref) > 0) {
        control_value <- found_pref[1]
      } else if (length(group_vals) == 2) {
        # pick the first sorted value for reproducibility
        control_value <- sort(group_vals)[1]
      } else {
        # fallback: pick the first value
        control_value <- group_vals[1]
      }
      if (args$verbose) cat("Auto-selected control group:", control_value, "\n", file = stderr())
    }

    diff_args <- list(
      x = result_se,
      samples = group_col,
      control = control_value,
      method = args$diff_method,
      test = args$test,
      randomizations = args$randomizations,
      pcorr = args$pcorr,
      verbose = args$verbose
    )

    diff_args <- diff_args[!sapply(diff_args, is.null)]

    diff_res <- do.call(SplicingFactory::calculate_difference, diff_args)

    tryCatch({
      write.table(
        diff_res,
        file = args$out_diff,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      if (args$verbose) cat("  Saved differential analysis to:", args$out_diff, "\n", file = stderr())
    }, error = function(e) {
      warning("Failed to write differential analysis output: ", conditionMessage(e))
    })
  }
}

# ============================================================================
# Output SummarizedExperiment (optional)
# ============================================================================

if (args$output_summarizedExperiment_RDS && !is.null(args$out_se)) {
  saveRDS(result_se, file = args$out_se)
  if (args$verbose) {
    cat(
      "  Saved SummarizedExperiment to:",
      args$out_se,
      "\n",
      file = stderr()
    )
  }
}

# ============================================================================
# Session info
# ============================================================================

if (!is.null(args$session_out)) {
  sink(file = args$session_out)
  print(sessionInfo())
  sink()
  if (args$verbose) {
    cat("  Saved session info to:", args$session_out, "\n", file = stderr())
  }
}

if (args$verbose) {
  cat("SplicingFactory analysis complete!\n", file = stderr())
}

q(status = 0)
