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
  "--assayno",
  help = "Assay number for SummarizedExperiment input (default: 1)",
  type = "integer",
  default = 1,
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
  "--group_column",
  help = "Column name in colData for grouping samples in boxplot",
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

if (grepl("\\.rds$", args$input, ignore.case = TRUE)) {
  # RDS input (SummarizedExperiment or tximport)
  input_data <- readRDS(args$input)

  if (is(input_data, "SummarizedExperiment")) {
    # Direct SummarizedExperiment
    se_data <- input_data
    readcounts <- assay(se_data, i = args$assayno)
    genes <- rownames(se_data)
  } else if (is.list(input_data)) {
    # tximport list
    if (use_tpm && !is.null(input_data$abundance)) {
      readcounts <- input_data$abundance
    } else if (!is.null(input_data$counts)) {
      readcounts <- input_data$counts
    } else {
      stop("RDS missing 'counts' or 'abundance' for tximport")
    }
    genes <- rownames(readcounts)
  } else {
    stop("RDS input must be either SummarizedExperiment or tximport list")
  }
} else {
  # TSV input
  if (is.null(args$genes)) {
    stop("--genes parameter required for TSV input")
  }

  # Read transcript counts matrix
  readcounts <- as.matrix(read.table(
    args$input,
    sep = "\t",
    header = TRUE,
    row.names = 1
  ))
  readcounts <- apply(readcounts, 2, as.numeric) # Ensure numeric type
  rownames(readcounts) <- rownames(as.matrix(read.table(
    args$input,
    sep = "\t",
    header = TRUE,
    row.names = 1
  )))

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
    # Single-column file: treat as gene vector aligned to rows of the matrix
    gene_vec <- as.character(gene_map_raw[[1]])
    if (length(gene_vec) != nrow(readcounts)) {
      msg <- paste(
        "Gene mapping file has 1 column but length does not match",
        "number of transcripts"
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

# ============================================================================
# Generate plots
# ============================================================================

if (!is.null(args$plot_dir) && dir.exists(args$plot_dir)) {
  if (args$verbose) {
    cat("Generating plots in format:",
      args$plot_format,
      "\n",
      file = stderr()
    )
  }

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
  plot_data_long <- pivot_longer(
    plot_data,
    -Gene,
    names_to = "Sample",
    values_to = "Diversity"
  )

  # Add grouping information if provided
  if (!is.null(args$group_column) && args$group_column != "") {
    cd <- colData(result_se)
    if (!is.null(cd) && args$group_column %in% colnames(cd)) {
      group_vec <- cd[[args$group_column]]
      plot_data_long$Group <- group_vec[
        match(plot_data_long$Sample, colnames(result_se))
      ]
    }
  }

  # Density plot
  p_density <- ggplot(plot_data_long, aes(x = Diversity, fill = Sample)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(
      title = paste("Diversity Distribution -", args$method),
      x = "Diversity Value",
      y = "Density"
    )

  # Boxplot
  has_group_plot <- (
    !is.null(args$group_column) && args$group_column != "" &&
      "Group" %in% colnames(plot_data_long)
  )

  if (has_group_plot) {
    p_boxplot <- ggplot(
      plot_data_long,
      aes(x = Group, y = Diversity, fill = Group)
    ) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      labs(
        title = paste("Diversity by", args$group_column),
        x = args$group_column,
        y = "Diversity Value"
      )
  } else {
    p_boxplot <- ggplot(
      plot_data_long,
      aes(x = Sample, y = Diversity, fill = Sample)
    ) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = paste("Diversity by Sample -", args$method),
        x = "Sample",
        y = "Diversity Value"
      )
  }

  # Heatmap using geom_tile
  plot_data_heatmap <- plot_data_long
  p_heatmap <- ggplot(
    plot_data_heatmap,
    aes(x = Sample, y = Gene, fill = Diversity)
  ) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste("Diversity Heatmap -", args$method),
      x = "Sample",
      y = "Gene"
    )

  # Save plots - ensure proper file naming based on format
  if (args$plot_format == "pdf") {
    ggsave(
      file.path(args$plot_dir, "diversity_density.pdf"),
      plot = p_density,
      device = "pdf",
      width = 8,
      height = 6
    )
    ggsave(
      file.path(args$plot_dir, "diversity_boxplot.pdf"),
      plot = p_boxplot,
      device = "pdf",
      width = 8,
      height = 6
    )
    ggsave(
      file.path(args$plot_dir, "diversity_heatmap.pdf"),
      plot = p_heatmap,
      device = "pdf",
      width = 10,
      height = 8
    )
  } else {
    ggsave(
      file.path(args$plot_dir, "diversity_density.png"),
      plot = p_density,
      device = "png",
      width = 8,
      height = 6,
      dpi = 100
    )
    ggsave(
      file.path(args$plot_dir, "diversity_boxplot.png"),
      plot = p_boxplot,
      device = "png",
      width = 8,
      height = 6,
      dpi = 100
    )
    ggsave(
      file.path(args$plot_dir, "diversity_heatmap.png"),
      plot = p_heatmap,
      device = "png",
      width = 10,
      height = 8,
      dpi = 100
    )
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

# ============================================================================
# Output SummarizedExperiment (optional)
# ============================================================================

if (save_se_requested && !is.null(args$out_se)) {
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
