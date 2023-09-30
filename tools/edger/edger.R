# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: htmlPath", "R", 1, "character"      -Path to html file linking to other outputs
#       outPath", "o", 1, "character"       -Path to folder to write all output to
#       filesPath", "j", 2, "character"     -JSON list object if multiple files input
#       matrixPath", "m", 2, "character"    -Path to count matrix
#       factFile", "f", 2, "character"      -Path to factor information file
#       factInput", "i", 2, "character"     -String containing factors if manually input
#       formula", "F", 2, "character".      -String containing a formula to override default use of factInput
#       annoPath", "a", 2, "character"      -Path to input containing gene annotations
#       contrastData", "C", 1, "character"  -String containing contrasts of interest
#       cpmReq", "c", 2, "double"           -Float specifying cpm requirement
#       cntReq", "z", 2, "integer"          -Integer specifying minimum total count requirement
#       sampleReq", "s", 2, "integer"       -Integer specifying cpm requirement
#       normCounts", "x", 0, "logical"      -String specifying if normalised counts should be output
#       rdaOpt", "r", 0, "logical"          -String specifying if RData should be output
#       lfcReq", "l", 1, "double"           -Float specifying the log-fold-change requirement
#       pValReq", "p", 1, "double"          -Float specifying the p-value requirement
#       pAdjOpt", "d", 1, "character"       -String specifying the p-value adjustment method
#       normOpt", "n", 1, "character"       -String specifying type of normalisation used
#       robOpt", "b", 0, "logical"          -String specifying if robust options should be used
#       lrtOpt", "t", 0, "logical"          -String specifying whether to perform LRT test instead
#
# OUT:
#       MDS Plot
#       BCV Plot
#       QL Plot
#       MD Plot
#       Expression Table
#       HTML file linking to the ouputs
# Optional:
#       Normalised counts Table
#       RData file
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Oct 2017 (some code taken from the DESeq2 wrapper)

# Record starting time
time_start <- as.character(Sys.time())

# setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, FALSE)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Load all required libraries
library(methods, quietly = TRUE, warn.conflicts = FALSE)
library(statmod, quietly = TRUE, warn.conflicts = FALSE)
library(splines, quietly = TRUE, warn.conflicts = FALSE)
library(edgeR, quietly = TRUE, warn.conflicts = FALSE)
library(limma, quietly = TRUE, warn.conflicts = FALSE)
library(scales, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)

################################################################################
### Function Delcaration
################################################################################
# Function to sanitise contrast equations so there are no whitespaces
# surrounding the arithmetic operators, leading or trailing whitespace
sanitise_equation <- function(equation) {
  equation <- gsub(" *[+] *", "+", equation)
  equation <- gsub(" *[-] *", "-", equation)
  equation <- gsub(" *[/] *", "/", equation)
  equation <- gsub(" *[*] *", "*", equation)
  equation <- gsub("^\\s+|\\s+$", "", equation)
  return(equation)
}

# Function to sanitise group information
sanitise_groups <- function(string) {
  string <- gsub(" *[,] *", ",", string)
  string <- gsub("^\\s+|\\s+$", "", string)
  return(string)
}

# Function to change periods to whitespace in a string
unmake_names <- function(string) {
  string <- gsub(".", " ", string, fixed = TRUE)
  return(string)
}

# Generate output folder and paths
make_out <- function(filename) {
  return(paste0(out_path, "/", filename))
}

# Generating design information
paste_listname <- function(string) {
  return(paste0("factors$", string))
}

# Create cata function: default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)
cata <- function(..., file = opt$htmlPath, sep = "", fill = FALSE, labels = NULL,
                 append = TRUE) {
  if (is.character(file)) {
    if (file == "") {
      file <- stdout()
    } else if (substring(file, 1L, 1L) == "|") {
      file <- pipe(substring(file, 2L), "w")
      on.exit(close(file))
    } else {
      file <- file(file, ifelse(append, "a", "w"))
      on.exit(close(file))
    }
  }
  .Internal(cat(list(...), file, sep, fill, labels, append))
}

# Function to write code for html head and title
html_head <- function(title) {
  cata("<head>\n")
  cata("<title>", title, "</title>\n")
  cata("</head>\n")
}

# Function to write code for html links
html_link <- function(address, label = address) {
  cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

# Function to write code for html images
html_image <- function(source, label = source, height = 600, width = 600) {
  cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
  cata("\" width=\"", width, "\"/>\n")
}

# Function to write code for html list items
list_item <- function(...) {
  cata("<li>", ..., "</li>\n")
}

table_item <- function(...) {
  cata("<td>", ..., "</td>\n")
}

table_head_item <- function(...) {
  cata("<th>", ..., "</th>\n")
}

################################################################################
### Input Processing
################################################################################

# Collect arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "htmlPath", "R", 1, "character",
  "outPath", "o", 1, "character",
  "filesPath", "j", 2, "character",
  "matrixPath", "m", 2, "character",
  "factFile", "f", 2, "character",
  "formula", "F", 2, "character",
  "factInput", "i", 2, "character",
  "annoPath", "a", 2, "character",
  "contrastData", "C", 1, "character",
  "cpmReq", "c", 1, "double",
  "totReq", "y", 0, "logical",
  "cntReq", "z", 1, "integer",
  "sampleReq", "s", 1, "integer",
  "normCounts", "x", 0, "logical",
  "rdaOpt", "r", 0, "logical",
  "lfcReq", "l", 1, "double",
  "pValReq", "p", 1, "double",
  "pAdjOpt", "d", 1, "character",
  "normOpt", "n", 1, "character",
  "robOpt", "b", 0, "logical",
  "lrtOpt", "t", 0, "logical"
),
byrow = TRUE, ncol = 4
)
opt <- getopt(spec)


if (is.null(opt$matrixPath) && is.null(opt$filesPath)) {
  cat("A counts matrix (or a set of counts files) is required.\n")
  q(status = 1)
}

if (is.null(opt$cpmReq)) {
  filt_cpm <- FALSE
} else {
  filt_cpm <- TRUE
}

if (is.null(opt$cntReq) || is.null(opt$sampleReq)) {
  filt_smpcount <- FALSE
} else {
  filt_smpcount <- TRUE
}

if (is.null(opt$totReq)) {
  filt_totcount <- FALSE
} else {
  filt_totcount <- TRUE
}

if (is.null(opt$lrtOpt)) {
  want_lrt <- FALSE
} else {
  want_lrt <- TRUE
}

if (is.null(opt$rdaOpt)) {
  want_rda <- FALSE
} else {
  want_rda <- TRUE
}

if (is.null(opt$annoPath)) {
  have_anno <- FALSE
} else {
  have_anno <- TRUE
}

if (is.null(opt$normCounts)) {
  want_norm <- FALSE
} else {
  want_norm <- TRUE
}

if (is.null(opt$robOpt)) {
  want_robust <- FALSE
} else {
  want_robust <- TRUE
}


if (!is.null(opt$filesPath)) {
  # Process the separate count files (adapted from DESeq2 wrapper)
  library("rjson")
  parser <- newJSONParser()
  parser$addData(opt$filesPath)
  factor_list <- parser$getObject()
  factors <- sapply(factor_list, function(x) x[[1]])
  filenames_in <- unname(unlist(factor_list[[1]][[2]]))
  sampletable <- data.frame(
    sample = basename(filenames_in),
    filename = filenames_in,
    row.names = filenames_in,
    stringsAsFactors = FALSE
  )
  for (factor in factor_list) {
    factorname  <- factor[[1]]
    sampletable[[factorname]] <- character(nrow(sampletable))
    lvls <- sapply(factor[[2]], function(x) names(x))
    for (i in seq_along(factor[[2]])) {
      files <- factor[[2]][[i]][[1]]
      sampletable[files, factorname] <- lvls[i]
    }
    sampletable[[factorname]] <- factor(sampletable[[factorname]], levels = lvls)
  }
  rownames(sampletable) <- sampletable$sample
  rem <- c("sample", "filename")
  factors <- sampletable[, !(names(sampletable) %in% rem), drop = FALSE]

  # read in count files and create single table
  countfiles <- lapply(sampletable$filename, function(x) {
    read.delim(x, row.names = 1)
  })
  counts <- do.call("cbind", countfiles)
} else {
  # Process the single count matrix
  counts <- read.table(opt$matrixPath, header = TRUE, sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE)
  row.names(counts) <- counts[, 1]
  counts <- counts[, -1]
  countsrows <- nrow(counts)

  # Process factors
  if (is.null(opt$factInput)) {
    factordata <- read.table(opt$factFile, header = TRUE, sep = "\t", strip.white = TRUE, stringsAsFactors = TRUE)
    # check samples names match
    if (!any(factordata[, 1] %in% colnames(counts))) {
      stop("Sample IDs in factors file and count matrix don't match")
    }
    # order samples as in counts matrix
    factordata <- factordata[match(colnames(counts), factordata[, 1]), ]
    factors <- data.frame(sapply(factordata[, -1, drop = FALSE], make.names))
  } else {
    factors <- unlist(strsplit(opt$factInput, "|", fixed = TRUE))
    factordata <- list()
    for (fact in factors) {
      newfact <- unlist(strsplit(fact, split = "::"))
      factordata <- rbind(factordata, newfact)
    } # Factors have the form: FACT_NAME::LEVEL,LEVEL,LEVEL,LEVEL,... The first factor is the Primary Factor.

    # Set the row names to be the name of the factor and delete first row
    row.names(factordata) <- factordata[, 1]
    factordata <- factordata[, -1]
    factordata <- sapply(factordata, sanitise_groups)
    factordata <- sapply(factordata, strsplit, split = ",")
    factordata <- sapply(factordata, make.names)
    # Transform factor data into data frame of R factor objects
    factors <- data.frame(factordata)
  }
}

# if annotation file provided
if (have_anno) {
  geneanno <- read.table(opt$annoPath, header = TRUE, sep = "\t", quote = "", strip.white = TRUE, stringsAsFactors = FALSE)
}

# Create output directory
out_path <- opt$outPath
dir.create(out_path, showWarnings = FALSE)

# Check if contrastData is a file or not
if (file.exists(opt$contrastData)) {
  contrast_data <- unlist(read.table(opt$contrastData, sep="\t", header=TRUE)[[1]])
} else {
  # Split up contrasts separated by comma into a vector then sanitise
  contrast_data <- unlist(strsplit(opt$contrastData, split = ","))
}
contrast_data <- sanitise_equation(contrast_data)
contrast_data <- gsub(" ", ".", contrast_data, fixed = TRUE)

bcv_pdf <- make_out("bcvplot.pdf")
bcv_png <- make_out("bcvplot.png")
ql_pdf <- make_out("qlplot.pdf")
ql_png <- make_out("qlplot.png")
mds_pdf <- character() # Initialise character vector
mds_png <- character()
for (i in seq_len(ncol(factors))) {
  mds_pdf[i] <- make_out(paste0("mdsplot_", names(factors)[i], ".pdf"))
  mds_png[i] <- make_out(paste0("mdsplot_", names(factors)[i], ".png"))
}
md_pdf <- character()
md_png <- character()
top_out <- character()
for (i in seq_along(contrast_data)) {
  md_pdf[i] <- make_out(paste0("mdplot_", contrast_data[i], ".pdf"))
  md_png[i] <- make_out(paste0("mdplot_", contrast_data[i], ".png"))
  top_out[i] <- make_out(paste0("edgeR_", contrast_data[i], ".tsv"))
} # Save output paths for each contrast as vectors
norm_out <- make_out("edgeR_normcounts.tsv")
rda_out <- make_out("edgeR_analysis.RData")
session_out <- make_out("session_info.txt")

# Initialise data for html links and images, data frame with columns Label and
# Link
link_data <- data.frame(Label = character(), Link = character(), stringsAsFactors = FALSE)
image_data <- data.frame(Label = character(), Link = character(), stringsAsFactors = FALSE)

# Initialise vectors for storage of up/down/neutral regulated counts
up_count <- numeric()
down_count <- numeric()
flat_count <- numeric()

################################################################################
### Data Processing
################################################################################

# Extract counts and annotation data
data <- list()
data$counts <- counts
if (have_anno) {
  # order annotation by genes in counts (assumes gene ids are in 1st column of geneanno)
  annoord <- geneanno[match(row.names(counts), geneanno[, 1]), ]
  data$genes <- annoord
} else {
  data$genes <- data.frame(GeneID = row.names(counts))
}

# If filter crieteria set, filter out genes that do not have a required cpm/counts in a required number of
# samples. Default is no filtering
prefilter_count <- nrow(data$counts)

if (filt_cpm || filt_smpcount || filt_totcount) {
  if (filt_totcount) {
    keep <- rowSums(data$counts) >= opt$cntReq
  } else if (filt_smpcount) {
    keep <- rowSums(data$counts >= opt$cntReq) >= opt$sampleReq
  } else if (filt_cpm) {
    keep <- rowSums(cpm(data$counts) >= opt$cpmReq) >= opt$sampleReq
  }

  data$counts <- data$counts[keep, ]
  data$genes <- data$genes[keep, , drop = FALSE]
}

postfilter_count <- nrow(data$counts)
filtered_count <- prefilter_count - postfilter_count

# Name rows of factors according to their sample
row.names(factors) <- names(data$counts)
factor_list <- names(factors)

# Generating the DGEList object "data"
samplenames <- colnames(data$counts)
genes <- data$genes
data <- DGEList(data$counts)
colnames(data) <- samplenames
data$samples <- factors
data$genes <- genes


if (!is.null(opt$formula)) {
  formula <- opt$formula
  # sanitisation can be getting rid of the "~"
  if(!startsWith(formula, "~")) {
    formula <- paste0("~", formula)
  }
} else {
  formula <- "~0"
  for (i in seq_along(factor_list)) {
    formula <- paste(formula, factor_list[i], sep = "+")
  }
}

formula <- formula(formula)
design <- model.matrix(formula, factors)

for (i in seq_along(factor_list)) {
  colnames(design) <- gsub(factor_list[i], "", colnames(design), fixed = TRUE)
}

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method = opt$normOpt)

if (want_robust) {
  data <- estimateDisp(data, design = design, robust = TRUE)
} else {
  data <- estimateDisp(data, design = design)
}

# Generate contrasts information
contrasts <- makeContrasts(contrasts = contrast_data, levels = design)

################################################################################
### Data Output
################################################################################

# Plot MDS
labels <- names(counts)

# MDS plot
png(mds_png, width = 600, height = 600)
plotMDS(data, labels = labels, col = as.numeric(factors[, 1]), cex = 0.8, main = paste("MDS Plot:", names(factors)[1]))
img_name <- paste0("MDS Plot_", names(factors)[1], ".png")
img_addr <- paste0("mdsplot_", names(factors)[1], ".png")
image_data[1, ] <- c(img_name, img_addr)
invisible(dev.off())

pdf(mds_pdf)
plotMDS(data, labels = labels, col = as.numeric(factors[, 1]), cex = 0.8, main = paste("MDS Plot:", names(factors)[1]))
link_name <- paste0("MDS Plot_", names(factors)[1], ".pdf")
link_addr <- paste0("mdsplot_", names(factors)[1], ".pdf")
link_data[1, ] <- c(link_name, link_addr)
invisible(dev.off())

# If additional factors create additional MDS plots coloured by factor
if (ncol(factors) > 1) {
  for (i in 2:ncol(factors)) {
    png(mds_png[i], width = 600, height = 600)
    plotMDS(data, labels = labels, col = as.numeric(factors[, i]), cex = 0.8, main = paste("MDS Plot:", names(factors)[i]))
    img_name <- paste0("MDS Plot_", names(factors)[i], ".png")
    img_addr <- paste0("mdsplot_", names(factors)[i], ".png")
    image_data <- rbind(image_data, c(img_name, img_addr))
    invisible(dev.off())

    pdf(mds_pdf[i])
    plotMDS(data, labels = labels, col = as.numeric(factors[, i]), cex = 0.8, main = paste("MDS Plot:", names(factors)[i]))
    link_name <- paste0("MDS Plot_", names(factors)[i], ".pdf")
    link_addr <- paste0("mdsplot_", names(factors)[i], ".pdf")
    link_data <- rbind(link_data, c(link_name, link_addr))
    invisible(dev.off())
  }
}

# BCV Plot
png(bcv_png, width = 600, height = 600)
plotBCV(data, main = "BCV Plot")
img_name <- "BCV Plot"
img_addr <- "bcvplot.png"
image_data <- rbind(image_data, c(img_name, img_addr))
invisible(dev.off())

pdf(bcv_pdf)
plotBCV(data, main = "BCV Plot")
link_name <- paste0("BCV Plot.pdf")
link_addr <- paste0("bcvplot.pdf")
link_data <- rbind(link_data, c(link_name, link_addr))
invisible(dev.off())

# Generate fit
if (want_lrt) {
  fit <- glmFit(data, design)
} else {
  if (want_robust) {
    fit <- glmQLFit(data, design, robust = TRUE)
  } else {
    fit <- glmQLFit(data, design)
  }

  # Plot QL dispersions
  png(ql_png, width = 600, height = 600)
  plotQLDisp(fit, main = "QL Plot")
  img_name <- "QL Plot"
  img_addr <- "qlplot.png"
  image_data <- rbind(image_data, c(img_name, img_addr))
  invisible(dev.off())

  pdf(ql_pdf)
  plotQLDisp(fit, main = "QL Plot")
  link_name <- "QL Plot.pdf"
  link_addr <- "qlplot.pdf"
  link_data <- rbind(link_data, c(link_name, link_addr))
  invisible(dev.off())
}

# Save normalised counts (log2cpm)
if (want_norm) {
  normalised_counts <- cpm(data, normalized.lib.sizes = TRUE, log = TRUE)
  normalised_counts <- data.frame(data$genes, normalised_counts)
  write.table(normalised_counts, file = norm_out, row.names = FALSE, sep = "\t", quote = FALSE)
  link_data <- rbind(link_data, c("edgeR_normcounts.tsv", "edgeR_normcounts.tsv"))
}


for (i in seq_along(contrast_data)) {
  if (want_lrt) {
    res <- glmLRT(fit, contrast = contrasts[, i])
  } else {
    res <- glmQLFTest(fit, contrast = contrasts[, i])
  }

  status <- decideTestsDGE(res,
    adjust.method = opt$pAdjOpt, p.value = opt$pValReq,
    lfc = opt$lfcReq
  )
  sum_status <- summary(status)

  # Collect counts for differential expression
  up_count[i] <- sum_status["Up", ]
  down_count[i] <- sum_status["Down", ]
  flat_count[i] <- sum_status["NotSig", ]

  # Write top expressions table
  top <- topTags(res, adjust.method = opt$pAdjOpt, n = Inf, sort.by = "PValue")
  write.table(top, file = top_out[i], row.names = FALSE, sep = "\t", quote = FALSE)

  link_name <- paste0("edgeR_", contrast_data[i], ".tsv")
  link_addr <- paste0("edgeR_", contrast_data[i], ".tsv")
  link_data <- rbind(link_data, c(link_name, link_addr))

  # Plot MD (log ratios vs mean difference) using limma package
  pdf(md_pdf[i])
  limma::plotMD(res,
    status = status,
    main = paste("MD Plot:", unmake_names(contrast_data[i])),
    hl.col = alpha(c("firebrick", "blue"), 0.4), values = c(1, -1),
    xlab = "Average Expression", ylab = "logFC"
  )

  abline(h = 0, col = "grey", lty = 2)

  link_name <- paste0("MD Plot_", contrast_data[i], ".pdf")
  link_addr <- paste0("mdplot_", contrast_data[i], ".pdf")
  link_data <- rbind(link_data, c(link_name, link_addr))
  invisible(dev.off())

  png(md_png[i], height = 600, width = 600)
  limma::plotMD(res,
    status = status,
    main = paste("MD Plot:", unmake_names(contrast_data[i])),
    hl.col = alpha(c("firebrick", "blue"), 0.4), values = c(1, -1),
    xlab = "Average Expression", ylab = "logFC"
  )

  abline(h = 0, col = "grey", lty = 2)

  img_name <- paste0("MD Plot_", contrast_data[i], ".png")
  img_addr <- paste0("mdplot_", contrast_data[i], ".png")
  image_data <- rbind(image_data, c(img_name, img_addr))
  invisible(dev.off())
}
sig_diff <- data.frame(Up = up_count, Flat = flat_count, Down = down_count)
row.names(sig_diff) <- contrast_data

# Save relevant items as rda object
if (want_rda) {
  if (want_norm) {
    save(counts, data, status, normalised_counts, labels, factors, fit, res, top, contrasts, design,
      file = rda_out, ascii = TRUE
    )
  } else {
    save(counts, data, status, labels, factors, fit, res, top, contrasts, design,
      file = rda_out, ascii = TRUE
    )
  }
  link_data <- rbind(link_data, c("edgeR_analysis.RData", "edgeR_analysis.RData"))
}

# Record session info
writeLines(capture.output(sessionInfo()), session_out)
link_data <- rbind(link_data, c("Session Info", "session_info.txt"))

# Record ending time and calculate total run time
time_end <- as.character(Sys.time())
time_taken <- capture.output(round(difftime(time_end, time_start), digits = 3))
time_taken <- gsub("Time difference of ", "", time_taken, fixed = TRUE)

################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file = opt$htmlPath)

cata("<html>\n")

cata("<body>\n")
cata("<h3>edgeR Analysis Output:</h3>\n")
cata("Links to PDF copies of plots are in 'Plots' section below.<br />\n")

html_image(image_data$Link[1], image_data$Label[1])

for (i in 2:nrow(image_data)) {
  html_image(image_data$Link[i], image_data$Label[i])
}

cata("<h4>Differential Expression Counts:</h4>\n")

cata("<table border=\"1\" cellpadding=\"4\">\n")
cata("<tr>\n")
table_item()
for (i in colnames(sig_diff)) {
  table_head_item(i)
}
cata("</tr>\n")
for (i in seq_len(nrow(sig_diff))) {
  cata("<tr>\n")
  table_head_item(unmake_names(row.names(sig_diff)[i]))
  for (j in seq_len(ncol(sig_diff))) {
    table_item(as.character(sig_diff[i, j]))
  }
  cata("</tr>\n")
}
cata("</table>")

cata("<h4>Plots:</h4>\n")
for (i in seq_len(nrow(link_data))) {
  if (grepl(".pdf", link_data$Link[i])) {
    html_link(link_data$Link[i], link_data$Label[i])
  }
}

cata("<h4>Tables:</h4>\n")
for (i in seq_len(nrow(link_data))) {
  if (grepl(".tsv", link_data$Link[i])) {
    html_link(link_data$Link[i], link_data$Label[i])
  }
}

if (want_rda) {
  cata("<h4>R Data Objects:</h4>\n")
  for (i in seq_len(nrow(link_data))) {
    if (grepl(".RData", link_data$Link[i])) {
      html_link(link_data$Link[i], link_data$Label[i])
    }
  }
}

cata("<p>Alt-click links to download file.</p>\n")
cata("<p>Click floppy disc icon associated history item to download ")
cata("all files.</p>\n")
cata("<p>.tsv files can be viewed in Excel or any spreadsheet program.</p>\n")

cata("<h4>Additional Information</h4>\n")
cata("<ul>\n")

if (filt_cpm || filt_smpcount || filt_totcount) {
  if (filt_cpm) {
    temp_str <- paste(
      "Genes without more than", opt$cpmReq,
      "CPM in at least", opt$sampleReq, "samples are insignificant",
      "and filtered out."
    )
  } else if (filt_smpcount) {
    temp_str <- paste(
      "Genes without more than", opt$cntReq,
      "counts in at least", opt$sampleReq, "samples are insignificant",
      "and filtered out."
    )
  } else if (filt_totcount) {
    temp_str <- paste(
      "Genes without more than", opt$cntReq,
      "counts, after summing counts for all samples, are insignificant",
      "and filtered out."
    )
  }

  list_item(temp_str)
  filter_prop <- round(filtered_count / prefilter_count * 100, digits = 2)
  temp_str <- paste0(
    filtered_count, " of ", prefilter_count, " (", filter_prop,
    "%) genes were filtered out for low expression."
  )
  list_item(temp_str)
}
list_item(opt$normOpt, " was the method used to normalise library sizes.")
if (want_lrt) {
  list_item("The edgeR likelihood ratio test was used.")
} else {
  if (want_robust) {
    list_item("The edgeR quasi-likelihood test was used with robust settings (robust=TRUE with estimateDisp and glmQLFit).")
  } else {
    list_item("The edgeR quasi-likelihood test was used.")
  }
}
if (opt$pAdjOpt != "none") {
  if (opt$pAdjOpt == "BH" || opt$pAdjOpt == "BY") {
    temp_str <- paste0(
      "MD-Plot highlighted genes are significant at FDR ",
      "of ", opt$pValReq, " and exhibit log2-fold-change of at ",
      "least ", opt$lfcReq, "."
    )
    list_item(temp_str)
  } else if (opt$pAdjOpt == "holm") {
    temp_str <- paste0(
      "MD-Plot highlighted genes are significant at adjusted ",
      "p-value of ", opt$pValReq, "  by the Holm(1979) ",
      "method, and exhibit log2-fold-change of at least ",
      opt$lfcReq, "."
    )
    list_item(temp_str)
  }
} else {
  temp_str <- paste0(
    "MD-Plot highlighted genes are significant at p-value ",
    "of ", opt$pValReq, " and exhibit log2-fold-change of at ",
    "least ", opt$lfcReq, "."
  )
  list_item(temp_str)
}
cata("</ul>\n")

cata("<h4>Summary of experimental data:</h4>\n")

cata("<p>*CHECK THAT SAMPLES ARE ASSOCIATED WITH CORRECT GROUP(S)*</p>\n")

cata("<table border=\"1\" cellpadding=\"3\">\n")
cata("<tr>\n")
table_head_item("SampleID")
table_head_item(names(factors)[1], " (Primary Factor)")

if (ncol(factors) > 1) {
  for (i in names(factors)[2:length(names(factors))]) {
    table_head_item(i)
  }
  cata("</tr>\n")
}

for (i in seq_len(nrow((factors)))) {
  cata("<tr>\n")
  table_head_item(row.names(factors)[i])
  for (j in seq_len(ncol(factors))) {
    table_item(as.character(unmake_names(factors[i, j])))
  }
  cata("</tr>\n")
}
cata("</table>")

for (i in seq_len(nrow(link_data))) {
  if (grepl("session_info", link_data$Link[i])) {
    html_link(link_data$Link[i], link_data$Label[i])
  }
}

cata("<table border=\"0\">\n")
cata("<tr>\n")
table_item("Task started at:")
table_item(time_start)
cata("</tr>\n")
cata("<tr>\n")
table_item("Task ended at:")
table_item(time_end)
cata("</tr>\n")
cata("<tr>\n")
table_item("Task run time:")
table_item(time_taken)
cata("<tr>\n")
cata("</table>\n")

cata("</body>\n")
cata("</html>")
