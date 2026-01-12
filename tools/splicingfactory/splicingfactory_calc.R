library("argparse", quietly = TRUE, warn.conflicts = FALSE)
library("SplicingFactory", quietly = TRUE, warn.conflicts = FALSE)
library("ggplot2", quietly = TRUE, warn.conflicts = FALSE)
library("SplicingFactory", quietly = TRUE, warn.conflicts = FALSE)
library("SummarizedExperiment", quietly = TRUE, warn.conflicts = FALSE)
library("reshape2", quietly = TRUE, warn.conflicts = FALSE)

options(
  show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  if (!suppressWarnings(requireNamespace("argparse", quietly = TRUE))) {
    cat("Required R package 'argparse' is not installed.\n", file = stderr())
    q(status = 1)
  }
  library(argparse)
})

parser <- ArgumentParser(
  description = "Improved Galaxy runner for SplicingFactory::calculate_diversity"
)

parser$add_argument(
  "-i", "--input",
  help = "Input file: TSV matrix or RDS (tximport list or SummarizedExperiment)",
  required = TRUE
)
parser$add_argument(
  "--input_type",
  help = "auto|matrix|tximport|rds (default: auto)",
  default = "auto"
)
parser$add_argument(
  "--tx2gene",
  help = paste(
    "Two-column TSV mapping transcripts to genes",
    "(transcript\tgene). Required for transcript-level matrix",
    "or tximport list unless RDS contains gene mapping"
  ),
  default = NULL
)
parser$add_argument(
  "--method",
  help = "Method: naive, laplace, gini, simpson, invsimpson",
  default = "laplace"
)
parser$add_argument(
  "-n", "--norm",
  help = "Normalize entropy: True/False",
  default = "True"
)
parser$add_argument(
  "-p", "--tpm",
  help = "Use TPM values from tximport list: True/False",
  default = "False"
)
parser$add_argument(
  "--assayno",
  help = "Assay number for SummarizedExperiment input (default: 1)",
  type = "integer",
  default = 1
)
parser$add_argument(
  "-o", "--out_div",
  help = "Output diversity TSV file",
  default = "diversity.tsv"
)
parser$add_argument(
  "-r", "--out_se",
  help = "Output SummarizedExperiment RDS file",
  default = "diversity_se.rds"
)
parser$add_argument(
  "--plot_dir",
  help = "Directory to save ggplot figures (if set). Files: diversity_density.png, diversity_boxplot.png, diversity_heatmap.png",
  default = NULL
)
parser$add_argument(
  "--plot_format",
  help = "Image format for plots: png|pdf (default: png)",
  default = "png"
)
parser$add_argument(
  "--group_column",
  help = "Column name in colData(res_se) to group samples for boxplot",
  default = NULL
)
parser$add_argument(
  "-s", "--session_out",
  help = "sessionInfo() output file",
  default = "session.txt"
)
parser$add_argument(
  "-v", "--verbose",
  help = "Verbose logging",
  action = "store_true",
  default = FALSE
)

# Diagnostic: print raw args (what R actually sees). This helps debug parse_args errors under planemo/galaxy.
cat("RAW_CMD_ARGS:", paste(shQuote(commandArgs(trailingOnly = TRUE)), collapse = " "), "\n", file = stderr())

args <- parser$parse_args()

verbose <- args$verbose
if (verbose) {
  cat("Input:", args$input, "type:", args$input_type, "\n")
}

# Detect the input type
input_type <- args$input_type
if (input_type == "auto") {
  # try to read RDS first
  rds_obj <- tryCatch(readRDS(args$input), error = function(e) NULL)
  if (!is.null(rds_obj)) {
    input_type <- "rds"
    if (verbose) {
      cat("Detected RDS input.\n")
    }
  } else {
    input_type <- "matrix"
    if (verbose) {
      cat("Assuming TSV matrix input.\n")
    }
  }
}

#  Check for tx2gene dependency based on input type
if (is.null(args$tx2gene)) {
  if (input_type == "matrix") {
    stop("For matrix input, --tx2gene is required")
  } else if (input_type == "tximport") {
    stop("tx2gene mapping is required for tximport list input")
  }
}

# helper to read TSV with optional gz
read_tsv_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path)) {
    gzfile(path, "rt")
  } else {
    path
  }
  df <- tryCatch(
    read.table(
      con,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      row.names = 1,
      stringsAsFactors = FALSE
    ),
    error = function(e) {
      stop("Failed to read TSV matrix: ", conditionMessage(e))
    }
  )
  as.matrix(df)
}

# helper to read tx2gene mapping
read_tx2gene <- function(path) {
  con <- if (grepl("\\.gz$", path)) {
    gzfile(path, "rt")
  } else {
    path
  }
  m <- tryCatch(
    read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) {
      stop("Failed to read tx2gene mapping: ", conditionMessage(e))
    }
  )
  if (ncol(m) < 2) {
    stop("tx2gene must have at least two columns: transcript and gene")
  }
  setNames(as.character(m[[2]]), as.character(m[[1]]))
}

expr <- NULL
genes <- NULL
tx2gene <- NULL

if (input_type == "matrix") {
  expr <- read_tsv_matrix(args$input)
  # tx2gene is required
  tx2gene <- read_tx2gene(args$tx2gene)
  tx_ids <- rownames(expr)
  if (is.null(tx_ids)) {
    stop("Input expression matrix must have transcript IDs as rownames")
  }
  genes <- tx2gene[tx_ids]
  if (any(is.na(genes))) {
    stop("Some transcripts in the expression matrix are missing in tx2gene mapping")
  }
} else if (input_type == "tximport") {
  # expect RDS containing tximport-style list
  obj <- tryCatch(
    readRDS(args$input),
    error = function(e) {
      stop("Failed to read RDS: ", conditionMessage(e))
    }
  )
  if (!is.list(obj) || !("counts" %in% names(obj))) {
    stop("RDS does not appear to contain a tximport list with 'counts'")
  }
  expr <- obj
  if (is.null(args$tx2gene)) {
    stop("tx2gene mapping is required for tximport list input")
  }
  tx2gene <- read_tx2gene(args$tx2gene)
  # will pass tximport list and supply genes vector later
} else if (input_type == "rds") {
  obj <- tryCatch(
    readRDS(args$input),
    error = function(e) {
      stop("Failed to read RDS: ", conditionMessage(e))
    }
  )
  # if SummarizedExperiment
  if (
    inherits(obj, "SummarizedExperiment") ||
      inherits(obj, "RangedSummarizedExperiment")
  ) {
    expr <- obj
  } else if (is.list(obj) && ("counts" %in% names(obj))) {
    expr <- obj
    if (!is.null(args$tx2gene)) {
      tx2gene <- read_tx2gene(args$tx2gene)
    }
  } else {
    stop(
      paste(
        "Unsupported RDS object. Expect a SummarizedExperiment",
        "or tximport list."
      )
    )
  }
} else {
  stop("Unsupported --input_type: ", input_type)
}

# convert norm/tpm strings to logical
norm_flag <- tolower(args$norm) %in% c("true", "t", "1")
tpm_flag <- tolower(args$tpm) %in% c("true", "t", "1")

if (verbose) {
  cat("Loading SplicingFactory and dependencies...\n")
}
pkg_load <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    cat(sprintf("Required R package '%s' is not installed.\n", pkg), file = stderr())
    q(status = 1)
  }
}

# Prepare genes vector for matrix/tximport cases
if (!is.null(expr) && is.matrix(expr)) {
  # already prepared above
  res_se <- tryCatch(
    {
      calculate_diversity(
        expr,
        genes = genes,
        method = args$method,
        norm = norm_flag,
        tpm = tpm_flag
      )
    },
    error = function(e) {
      stop(
        "Error running calculate_diversity on matrix: ",
        conditionMessage(e)
      )
    }
  )
} else if (!is.null(expr) && is.list(expr) && ("counts" %in% names(expr))) {
  # tximport list
  # ensure tx2gene exists
  if (is.null(args$tx2gene) && is.null(tx2gene)) {
    stop("tx2gene mapping is required for tximport inputs")
  }
  # calculate_diversity expects a genes vector matching transcripts; create it
  tx_ids <- rownames(expr$counts)
  if (is.null(tx_ids)) {
    stop("tximport counts must have rownames with transcript IDs")
  }
  if (is.null(tx2gene)) {
    tx2gene <- read_tx2gene(args$tx2gene)
  }
  genes <- tx2gene[tx_ids]
  if (any(is.na(genes))) {
    stop("Some transcripts in tximport counts are missing in tx2gene mapping")
  }
  res_se <- tryCatch(
    {
      calculate_diversity(
        expr,
        genes = genes,
        method = args$method,
        norm = norm_flag,
        tpm = tpm_flag
      )
    },
    error = function(e) {
      stop(
        "Error running calculate_diversity on tximport list: ",
        conditionMessage(e)
      )
    }
  )
} else if (
  !is.null(expr) && (
    inherits(expr, "SummarizedExperiment") ||
      inherits(expr, "RangedSummarizedExperiment")
  )
) {
  # SummarizedExperiment: pass directly, allow assay selection
  if (args$assayno <= 0) {
    stop("--assayno must be >= 1")
  }
  # calculate_diversity will extract the assay by assayno internally
  res_se <- tryCatch(
    {
      calculate_diversity(
        expr,
        genes = NULL,
        method = args$method,
        norm = norm_flag,
        tpm = tpm_flag,
        assayno = args$assayno
      )
    },
    error = function(e) {
      stop(
        "Error running calculate_diversity on SummarizedExperiment: ",
        conditionMessage(e)
      )
    }
  )
} else {
  stop("Unsupported input object type for calculation")
}

# write table: assays(diversity)
if (!("diversity" %in% SummarizedExperiment::assayNames(res_se))) {
  stop("Resulting SummarizedExperiment does not contain an assay named 'diversity'.")
}
div_mat <- SummarizedExperiment::assay(res_se, "diversity")
out_df <- data.frame(
  Gene = rownames(div_mat),
  div_mat,
  check.names = FALSE
)
write.table(
  out_df,
  file = args$out_div,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# save SummarizedExperiment as RDS
saveRDS(res_se, file = args$out_se)

# optional plotting similar to vignette
if (!is.null(args$plot_dir)) {
  if (!dir.exists(args$plot_dir)) dir.create(args$plot_dir, recursive = TRUE)
  fmt <- tolower(args$plot_format)
  open_dev <- function(path) {
    if (fmt == "pdf") pdf(path, width = 7, height = 5) else png(path, width = 900, height = 650, res = 120)
  }
  close_dev <- function() dev.off()

  div_df <- as.data.frame(div_mat)
  div_df$Gene <- rownames(div_mat)
  long_df <- reshape2::melt(div_df, id.vars = "Gene", variable.name = "Sample", value.name = "Diversity")

  # density plot across samples
  p1 <- ggplot(long_df, aes(x = Diversity, group = Sample)) +
    geom_density(alpha = 0.2, color = NA, fill = "steelblue") +
    labs(title = paste("Diversity (", args$method, ") density", sep = "")) +
    theme_minimal()
  open_dev(file.path(args$plot_dir, paste0("diversity_density.", fmt)))
  print(p1)
  close_dev()

  # boxplot grouped by condition if provided
  group_col <- NULL
  if (!is.null(args$group_column)) {
    cd <- tryCatch(as.data.frame(SummarizedExperiment::colData(res_se)), error = function(e) NULL)
    if (!is.null(cd) && args$group_column %in% colnames(cd)) {
      long_df$Group <- cd[match(long_df$Sample, rownames(cd)), args$group_column]
      group_col <- "Group"
    }
  }
  if (is.null(group_col)) {
    long_df$Group <- "all"
  }
  p2 <- ggplot(long_df, aes(x = Group, y = Diversity)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(title = paste("Diversity (", args$method, ") by group", sep = ""), x = ifelse(is.null(args$group_column), "All samples", args$group_column)) +
    theme_minimal()
  open_dev(file.path(args$plot_dir, paste0("diversity_boxplot.", fmt)))
  print(p2)
  close_dev()

  # heatmap of genes by samples
  mat <- as.matrix(div_mat)
  # cap extreme values for visualization
  finite_vals <- mat[is.finite(mat)]
  if (length(finite_vals) > 0) {
    q <- quantile(finite_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    mat[mat < q[1]] <- q[1]
    mat[mat > q[2]] <- q[2]
  }
  df_melt <- reshape2::melt(mat, varnames = c("Gene", "Sample"), value.name = "Diversity")
  p3 <- ggplot(df_melt, aes(Sample, Gene, fill = Diversity)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(title = paste("Diversity (", args$method, ") heatmap", sep = ""))
  open_dev(file.path(args$plot_dir, paste0("diversity_heatmap.", fmt)))
  print(p3)
  close_dev()
}

# session info
con <- file(args$session_out, open = "wt")
cat("SplicingFactory wrapper run\n", file = con)
cat("Arguments:\n", file = con)
arg_list <- as.list(args)
cat(
  paste(names(arg_list), unlist(arg_list), sep = "=", collapse = ", "),
  "\n",
  file = con
)
cat("\nSession info:\n", file = con)
capture.output(sessionInfo(), file = con)
close(con)

invisible(NULL)