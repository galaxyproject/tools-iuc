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
#       annoPath", "a", 2, "character"      -Path to input containing gene annotations
#       contrastFile", "C", 1, "character"  -Path to contrasts information file
#       contrastInput", "D", 1, "character" -String containing contrasts of interest
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
#       trend", "t", 1, "double"            -Float for prior.count if limma-trend is used instead of voom
#       weightOpt", "w", 0, "logical"       -String specifying if voomWithQualityWeights should be used
#       topgenes", "G", 1, "integer"        -Integer specifying no. of genes to highlight in volcano and heatmap
#       treatOpt", "T", 0, "logical"        -String specifying if TREAT function should be used
#       plots, "P", 1, "character"          -String specifying additional plots to be created
#
# OUT:
#       Density Plots (if filtering)
#       Box Plots (if normalising)
#       MDS Plot
#       Voom/SA plot
#       MD Plot
#       Volcano Plot
#       Heatmap
#       Expression Table
#       HTML file linking to the ouputs
# Optional:
#       Normalised counts Table
#       RData file
#
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Jun 2017, Jan 2018, May 2018

# Record starting time
time_start <- as.character(Sys.time())

# Load all required libraries
library(methods, quietly = TRUE, warn.conflicts = FALSE)
library(statmod, quietly = TRUE, warn.conflicts = FALSE)
library(splines, quietly = TRUE, warn.conflicts = FALSE)
library(edgeR, quietly = TRUE, warn.conflicts = FALSE)
library(limma, quietly = TRUE, warn.conflicts = FALSE)
library(scales, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)
library(gplots, quietly = TRUE, warn.conflicts = FALSE)

################################################################################
### Function Declaration
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

# Function to make contrast contain valid R names
sanitise_contrast <- function(string) {
    string <- strsplit(string, split = "-")
    string <- lapply(string, make.names)
    string <- lapply(string, paste, collapse = "-")
    return(string)
}

# Function to change periods to whitespace in a string
unmake_names <- function(string) {
    string <- gsub(".", " ", string, fixed = TRUE)
    return(string)
}

# Generate output folder and paths
make_out <- function(filename) {
    return(paste0(opt$outPath, "/", filename))
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
html_image <- function(source, label = source, height = 500, width = 500) {
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
spec <- matrix(
    c(
        "htmlPath", "R", 1, "character",
        "outPath", "o", 1, "character",
        "filesPath", "j", 2, "character",
        "matrixPath", "m", 2, "character",
        "factFile", "f", 2, "character",
        "factInput", "i", 2, "character",
        "annoPath", "a", 2, "character",
        "contrastFile", "C", 1, "character",
        "contrastInput", "D", 1, "character",
        "cpmReq", "c", 1, "double",
        "totReq", "y", 0, "logical",
        "cntReq", "z", 1, "integer",
        "sampleReq", "s", 1, "integer",
        "filtCounts", "F", 0, "logical",
        "normCounts", "x", 0, "logical",
        "rdaOpt", "r", 0, "logical",
        "lfcReq", "l", 1, "double",
        "pValReq", "p", 1, "double",
        "pAdjOpt", "d", 1, "character",
        "normOpt", "n", 1, "character",
        "robOpt", "b", 0, "logical",
        "trend", "t", 1, "double",
        "weightOpt", "w", 0, "logical",
        "topgenes", "G", 1, "integer",
        "treatOpt", "T", 0, "logical",
        "plots", "P", 1, "character",
        "libinfoOpt", "L", 0, "logical"
    ),
    byrow = TRUE, ncol = 4
)
opt <- getopt(spec)


if (is.null(opt$matrixPath) & is.null(opt$filesPath)) {
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

if (is.null(opt$filtCounts)) {
    want_filt <- FALSE
} else {
    want_filt <- TRUE
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

if (is.null(opt$weightOpt)) {
    want_weight <- FALSE
} else {
    want_weight <- TRUE
}

if (is.null(opt$trend)) {
    want_trend <- FALSE
    de_method <- "limma-voom"
} else {
    want_trend <- TRUE
    de_method <- "limma-trend"
    prior_count <- opt$trend
}

if (is.null(opt$treatOpt)) {
    want_treat <- FALSE
} else {
    want_treat <- TRUE
}

if (is.null(opt$libinfoOpt)) {
    want_libinfo <- FALSE
} else {
    want_libinfo <- TRUE
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
        factorname <- factor[[1]]
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
    counts <- read.table(opt$matrixPath, header = TRUE, sep = "\t", strip.white = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[, -1]
    countsrows <- nrow(counts)

    # Process factors
    if (is.null(opt$factInput)) {
        factordata <- read.table(opt$factFile, header = TRUE, sep = "\t", strip.white = TRUE, stringsAsFactors = TRUE)
        if (!setequal(factordata[, 1], colnames(counts))) {
            stop("Sample IDs in counts and factors files don't match")
        }
        # order samples as in counts matrix
        factordata <- factordata[match(colnames(counts), factordata[, 1]), ]
        factors <- factordata[, -1, drop = FALSE]
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
        # Transform factor data into data frame of R factor objects
        factors <- data.frame(factordata)
    }
}
# check there are the same number of samples in counts and factors
if (nrow(factors) != ncol(counts)) {
    stop("There are a different number of samples in the counts files and factors")
}
# make groups valid R names, required for makeContrasts
factors <- sapply(factors, make.names)
factors <- data.frame(factors, stringsAsFactors = TRUE)

# if annotation file provided
if (have_anno) {
    geneanno <- read.table(opt$annoPath, header = TRUE, sep = "\t", quote = "", strip.white = TRUE, stringsAsFactors = FALSE)
}

# Create output directory
dir.create(opt$outPath, showWarnings = FALSE)

# Process contrasts
if (is.null(opt$contrastInput)) {
    contrast_data <- read.table(opt$contrastFile, header = TRUE, sep = "\t", quote = "", strip.white = TRUE, stringsAsFactors = FALSE)
    contrast_data <- contrast_data[, 1, drop = TRUE]
} else {
    # Split up contrasts seperated by comma into a vector then sanitise
    contrast_data <- unlist(strsplit(opt$contrastInput, split = ","))
}
contrast_data <- sanitise_equation(contrast_data)
contrast_data <- gsub(" ", ".", contrast_data, fixed = TRUE)

# in case input groups start with numbers make the names valid R names, required for makeContrasts
cons <- NULL
cons_d <- NULL
for (i in contrast_data) {
    # if the contrast is a difference of differences e.g. (A-B)-(X-Y)
    if (grepl("\\)-\\(", i)) {
        i <- unlist(strsplit(i, split = "\\)-\\("))
        i <- gsub("\\(|\\)", "", i)
        for (j in i) {
            j <- sanitise_contrast(j)
            j <- paste0("(", j, ")")
            cons_d <- append(cons_d, unlist(j))
        }
        cons_d <- paste(cons_d, collapse = "-")
        cons <- append(cons, unlist(cons_d))
    } else {
        i <- sanitise_contrast(i)
        cons <- append(cons, unlist(i))
    }
}

plots <- character()
if (!is.null(opt$plots)) {
    plots <- unlist(strsplit(opt$plots, split = ","))
}

den_png <- make_out("densityplots.png")
den_pdf <- make_out("densityplots.pdf")
cpm_pdf <- make_out("cpmplots.pdf")
box_png <- make_out("boxplots.png")
box_pdf <- make_out("boxplots.pdf")
mdsscree_png <- make_out("mdsscree.png")
mdsscree_pdf <- make_out("mdsscree.pdf")
mdsx_pdf <- make_out("mdsplot_extra.pdf")
mdsx_png <- make_out("mdsplot_extra.png")
mdsam_pdf <- make_out("mdplots_samples.pdf")
md_pdf <- character() # Initialise character vector
vol_pdf <- character()
heat_pdf <- character()
strip_pdf <- character()
mdvol_png <- character()
top_out <- character()
glimma_out <- character()
for (i in seq_along(cons)) {
    con <- cons[i]
    con <- gsub("\\(|\\)", "", con)
    md_pdf[i] <- make_out(paste0("mdplot_", con, ".pdf"))
    vol_pdf[i] <- make_out(paste0("volplot_", con, ".pdf"))
    heat_pdf[i] <- make_out(paste0("heatmap_", con, ".pdf"))
    strip_pdf[i] <- make_out(paste0("stripcharts_", con, ".pdf"))
    mdvol_png[i] <- make_out(paste0("mdvolplot_", con, ".png"))
    top_out[i] <- make_out(paste0(de_method, "_", con, ".tsv"))
    glimma_out[i] <- make_out(paste0("glimma_", con, "/MD-Plot.html"))
}
filt_out <- make_out(paste0(de_method, "_", "filtcounts"))
norm_out <- make_out(paste0(de_method, "_", "normcounts"))
rda_out <- make_out(paste0(de_method, "_analysis.RData"))
session_out <- make_out("session_info.txt")

# Initialise data for html links and images, data frame with columns Label and
# Link
link_data <- data.frame(
    Label = character(), Link = character(),
    stringsAsFactors = FALSE
)
image_data <- data.frame(
    Label = character(), Link = character(),
    stringsAsFactors = FALSE
)

# Initialise vectors for storage of up/down/neutral regulated counts
up_count <- numeric()
down_count <- numeric()
flat_count <- numeric()

################################################################################
### Data Processing
################################################################################

# Extract counts and annotation data
print("Extracting counts")
data <- list()
data$counts <- counts
if (have_anno) {
    # order annotation by genes in counts (assumes gene ids are in 1st column of geneanno)
    annoord <- geneanno[match(row.names(counts), geneanno[, 1]), ]
    data$genes <- annoord
} else {
    data$genes <- data.frame(GeneID = row.names(counts))
}

# Creating naming data
samplenames <- colnames(data$counts)
sampleanno <- data.frame("sampleID" = samplenames, factors)
row.names(factors) <- samplenames # for "Summary of experimental data" table

# Creating colours for the groups
cols <- as.numeric(factors[, 1])
col.group <- palette()[cols]

# If filter crieteria set, filter out genes that do not have a required cpm/counts in a required number of
# samples. Default is no filtering
prefilter_count <- nrow(data$counts)
nsamples <- ncol(data$counts)

if (filt_cpm || filt_smpcount || filt_totcount) {
    if (filt_totcount) {
        keep <- rowSums(data$counts) >= opt$cntReq
    } else if (filt_smpcount) {
        keep <- rowSums(data$counts >= opt$cntReq) >= opt$sampleReq
    } else if (filt_cpm) {
        cpms <- cpm(data$counts)
        thresh <- cpms >= opt$cpmReq
        keep <- rowSums(thresh) >= opt$sampleReq

        if ("c" %in% plots) {
            # Plot CPM vs raw counts (to check threshold)
            pdf(cpm_pdf, width = 6.5, height = 10)
            par(mfrow = c(3, 2))
            for (i in seq_len(nsamples)) {
                plot(data$counts[, i], cpms[, i], xlim = c(0, 50), ylim = c(0, 3), main = samplenames[i], xlab = "Raw counts", ylab = "CPM")
                abline(v = 10, col = "red", lty = 2, lwd = 2)
                abline(h = opt$cpmReq, col = 4)
            }
            link_name <- "CpmPlots.pdf"
            link_addr <- "cpmplots.pdf"
            link_data <- rbind(link_data, data.frame(Label = link_name, Link = link_addr, stringsAsFactors = FALSE))
            invisible(dev.off())
        }
    }

    data$counts <- data$counts[keep, ]
    data$genes <- data$genes[keep, , drop = FALSE]

    if (want_filt) {
        print("Outputting filtered counts")
        filt_counts <- data.frame(data$genes, data$counts, check.names = FALSE)
        write.table(filt_counts, file = filt_out, row.names = FALSE, sep = "\t", quote = FALSE)
        link_data <- rbind(link_data, data.frame(Label = paste0(de_method, "_", "filtcounts.tsv"), Link = paste0(de_method, "_", "filtcounts"), stringsAsFactors = FALSE))
    }

    # Plot Density
    if ("d" %in% plots) {
        # PNG
        png(den_png, width = 1000, height = 500)
        par(mfrow = c(1, 2), cex.axis = 0.8)

        # before filtering
        lcpm1 <- cpm(counts, log = TRUE)
        plot(density(lcpm1[, 1]), col = col.group[1], lwd = 2, las = 2, main = "", xlab = "")
        title(main = "Density Plot: Raw counts", xlab = "Log-cpm")
        for (i in 2:nsamples) {
            den <- density(lcpm1[, i])
            lines(den$x, den$y, col = col.group[i], lwd = 2)
        }

        # after filtering
        lcpm2 <- cpm(data$counts, log = TRUE)
        plot(density(lcpm2[, 1]), col = col.group[1], lwd = 2, las = 2, main = "", xlab = "")
        title(main = "Density Plot: Filtered counts", xlab = "Log-cpm")
        for (i in 2:nsamples) {
            den <- density(lcpm2[, i])
            lines(den$x, den$y, col = col.group[i], lwd = 2)
        }
        legend("topright", samplenames, text.col = col.group, bty = "n")
        img_name <- "Densityplots.png"
        img_addr <- "densityplots.png"
        image_data <- rbind(image_data, data.frame(Label = img_name, Link = img_addr, stringsAsFactors = FALSE))
        invisible(dev.off())

        # PDF
        pdf(den_pdf, width = 14)
        par(mfrow = c(1, 2), cex.axis = 0.8)
        plot(density(lcpm1[, 1]), col = col.group[1], lwd = 2, las = 2, main = "", xlab = "")
        title(main = "Density Plot: Raw counts", xlab = "Log-cpm")
        for (i in 2:nsamples) {
            den <- density(lcpm1[, i])
            lines(den$x, den$y, col = col.group[i], lwd = 2)
        }
        plot(density(lcpm2[, 1]), col = col.group[1], lwd = 2, las = 2, main = "", xlab = "")
        title(main = "Density Plot: Filtered counts", xlab = "Log-cpm")
        for (i in 2:nsamples) {
            den <- density(lcpm2[, i])
            lines(den$x, den$y, col = col.group[i], lwd = 2)
        }
        legend("topright", samplenames, text.col = col.group, bty = "n")
        link_name <- "DensityPlots.pdf"
        link_addr <- "densityplots.pdf"
        link_data <- rbind(link_data, data.frame(Label = link_name, Link = link_addr, stringsAsFactors = FALSE))
        invisible(dev.off())
    }
}

postfilter_count <- nrow(data$counts)
filtered_count <- prefilter_count - postfilter_count

# Generating the DGEList object "y"
print("Generating DGEList object")
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts) # nolint
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
y <- new("DGEList", data)

print("Generating Design")
factor_list <- sapply(names(factors), paste_listname)
formula <- "~0"
for (i in seq_along(factor_list)) {
    formula <- paste(formula, factor_list[i], sep = "+")
}
formula <- formula(formula)
design <- model.matrix(formula)
for (i in seq_along(factor_list)) {
    colnames(design) <- gsub(factor_list[i], "", colnames(design), fixed = TRUE)
}

# Calculating normalising factors
print("Calculating Normalisation Factors")
logcounts <- y # store for plots
y <- calcNormFactors(y, method = opt$normOpt)

# Generate contrasts information
print("Generating Contrasts")
contrasts <- makeContrasts(contrasts = cons, levels = design)

################################################################################
### Data Output
################################################################################

# Plot Box plots (before and after normalisation)
if (opt$normOpt != "none" & "b" %in% plots) {
    png(box_png, width = 1000, height = 500)
    par(mfrow = c(1, 2), mar = c(6, 4, 2, 2) + 0.1)
    labels <- colnames(counts)

    lcpm1 <- cpm(y$counts, log = TRUE)
    boxplot(lcpm1, las = 2, col = col.group, xaxt = "n", xlab = "")
    axis(1, at = seq_along(labels), labels = FALSE)
    abline(h = median(lcpm1), col = 4)
    text(x = seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1, labels = labels, xpd = TRUE)
    title(main = "Box Plot: Unnormalised counts", ylab = "Log-cpm")

    lcpm2 <- cpm(y, log = TRUE)
    boxplot(lcpm2, las = 2, col = col.group, xaxt = "n", xlab = "")
    axis(1, at = seq_along(labels), labels = FALSE)
    text(x = seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1, labels = labels, xpd = TRUE)
    abline(h = median(lcpm2), col = 4)
    title(main = "Box Plot: Normalised counts", ylab = "Log-cpm")

    img_name <- "Boxplots.png"
    img_addr <- "boxplots.png"
    image_data <- rbind(image_data, data.frame(Label = img_name, Link = img_addr, stringsAsFactors = FALSE))
    invisible(dev.off())

    pdf(box_pdf, width = 14)
    par(mfrow = c(1, 2), mar = c(6, 4, 2, 2) + 0.1)
    boxplot(lcpm1, las = 2, col = col.group, xaxt = "n", xlab = "")
    axis(1, at = seq_along(labels), labels = FALSE)
    abline(h = median(lcpm1), col = 4)
    text(x = seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1, labels = labels, xpd = TRUE)
    title(main = "Box Plot: Unnormalised counts", ylab = "Log-cpm")
    boxplot(lcpm2, las = 2, col = col.group, xaxt = "n", xlab = "")
    axis(1, at = seq_along(labels), labels = FALSE)
    text(x = seq_along(labels), y = par("usr")[3] - 1, srt = 45, adj = 1, labels = labels, xpd = TRUE)
    abline(h = median(lcpm2), col = 4)
    title(main = "Box Plot: Normalised counts", ylab = "Log-cpm")
    link_name <- "BoxPlots.pdf"
    link_addr <- "boxplots.pdf"
    link_data <- rbind(link_data, data.frame(Label = link_name, Link = link_addr, stringsAsFactors = FALSE))
    invisible(dev.off())
}

# Plot MDS
print("Generating MDS plot")
labels <- names(counts)

# Scree plot (Variance Explained) code copied from Glimma

# get column of matrix
get_cols <- function(x, inds) {
    x[, inds, drop = FALSE]
}

x <- cpm(y, log = TRUE)
ndim <- nsamples - 1
nprobes <- nrow(x)
top <- 500
top <- min(top, nprobes)
cn <- colnames(x)
bad <- rowSums(is.finite(x)) < nsamples

if (any(bad)) {
    warning("Rows containing infinite values have been removed")
    x <- x[!bad, , drop = FALSE]
}

dd <- matrix(0, nrow = nsamples, ncol = nsamples, dimnames = list(cn, cn))
topindex <- nprobes - top + 1L
for (i in 2L:(nsamples)) {
    for (j in 1L:(i - 1L)) {
        dists <- (get_cols(x, i) - get_cols(x, j))^2
        dists <- sort.int(dists, partial = topindex)
        topdist <- dists[topindex:nprobes]
        dd[i, j] <- sqrt(mean(topdist))
    }
}

a1 <- suppressWarnings(cmdscale(as.dist(dd), k = min(ndim, 8), eig = TRUE))
eigen <- data.frame(name = 1:min(ndim, 8), eigen = round(a1$eig[1:min(ndim, 8)] / sum(a1$eig), 2))
png(mdsscree_png, width = 1000, height = 500)
par(mfrow = c(1, 2))
plotMDS(y, labels = samplenames, col = as.numeric(factors[, 1]), main = "MDS Plot: Dims 1 and 2")
barplot(eigen$eigen, names.arg = eigen$name, main = "Scree Plot: Variance Explained", xlab = "Dimension", ylab = "Proportion", las = 1)
img_name <- paste0("MDSPlot_", names(factors)[1], ".png")
img_addr <- "mdsscree.png"
image_data <- rbind(image_data, data.frame(Label = img_name, Link = img_addr, stringsAsFactors = FALSE))
invisible(dev.off())

pdf(mdsscree_pdf, width = 14)
par(mfrow = c(1, 2))
plotMDS(y, labels = samplenames, col = as.numeric(factors[, 1]), main = "MDS Plot: Dims 1 and 2")
barplot(eigen$eigen, names.arg = eigen$name, main = "Scree Plot: Variance Explained", xlab = "Dimension", ylab = "Proportion", las = 1)
link_name <- paste0("MDSPlot_", names(factors)[1], ".pdf")
link_addr <- "mdsscree.pdf"
link_data <- rbind(link_data, data.frame(Label = link_name, Link = link_addr, stringsAsFactors = FALSE))
invisible(dev.off())

# generate Glimma interactive MDS Plot
if ("i" %in% plots) {
    Glimma::glMDSPlot(y,
        labels = samplenames, groups = factors[, 1],
        folder = "glimma_MDS", launch = FALSE
    )
    link_name <- "Glimma_MDSPlot.html"
    link_addr <- "glimma_MDS/MDS-Plot.html"
    link_data <- rbind(link_data, c(link_name, link_addr))
}

if ("x" %in% plots) {
    png(mdsx_png, width = 1000, height = 500)
    par(mfrow = c(1, 2))
    for (i in 2:3) {
        dim1 <- i
        dim2 <- i + 1
        plotMDS(y, dim = c(dim1, dim2), labels = samplenames, col = as.numeric(factors[, 1]), main = paste("MDS Plot: Dims", dim1, "and", dim2))
    }
    img_name <- paste0("MDSPlot_extra.png")
    img_addr <- paste0("mdsplot_extra.png")
    image_data <- rbind(image_data, data.frame(Label = img_name, Link = img_addr, stringsAsFactors = FALSE))
    invisible(dev.off())

    pdf(mdsx_pdf, width = 14)
    par(mfrow = c(1, 2))
    for (i in 2:3) {
        dim1 <- i
        dim2 <- i + 1
        plotMDS(y, dim = c(dim1, dim2), labels = samplenames, col = as.numeric(factors[, 1]), main = paste("MDS Plot: Dims", dim1, "and", dim2))
    }
    link_name <- "MDSPlot_extra.pdf"
    link_addr <- "mdsplot_extra.pdf"
    link_data <- rbind(link_data, data.frame(Label = link_name, Link = link_addr, stringsAsFactors = FALSE))
    invisible(dev.off())
}

if ("m" %in% plots) {
    # Plot MD plots for individual samples
    print("Generating MD plots for samples")
    pdf(mdsam_pdf, width = 6.5, height = 10)
    par(mfrow = c(3, 2))
    for (i in 1:nsamples) {
        if (opt$normOpt != "none") {
            plotMD(logcounts, column = i, main = paste(colnames(logcounts)[i], "(before)"))
            abline(h = 0, col = "red", lty = 2, lwd = 2)
        }
        plotMD(y, column = i)
        abline(h = 0, col = "red", lty = 2, lwd = 2)
    }
    link_name <- "MDPlots_Samples.pdf"
    link_addr <- "mdplots_samples.pdf"
    link_data <- rbind(link_data, c(link_name, link_addr))
    invisible(dev.off())
}


if (want_trend) {
    # limma-trend approach
    logcpm <- cpm(y, log = TRUE, prior.count = opt$trend)
    fit <- lmFit(logcpm, design)
    fit$genes <- y$genes
    fit <- contrasts.fit(fit, contrasts)
    if (want_robust) {
        fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    } else {
        fit <- eBayes(fit, trend = TRUE, robust = FALSE)
    }

    plot_data <- logcpm

    # Save normalised counts (log2cpm)
    if (want_norm) {
        write.table(logcpm, file = norm_out, row.names = TRUE, sep = "\t", quote = FALSE)
        link_data <- rbind(link_data, c((paste0(de_method, "_", "normcounts.tsv")), (paste0(de_method, "_", "normcounts"))))
    }
} else {
    # limma-voom approach

    if (want_weight) {
        voomwts_pdf <- make_out("voomwtsplot.pdf")
        voomwts_png <- make_out("voomwtsplot.png")
        # Creating voom data object and plot
        png(voomwts_png, width = 1000, height = 500)
        vdata <- voomWithQualityWeights(y, design = design, plot = TRUE)
        img_name <- "VoomWithQualityWeightsPlot.png"
        img_addr <- "voomwtsplot.png"
        image_data <- rbind(image_data, c(img_name, img_addr))
        invisible(dev.off())

        pdf(voomwts_pdf, width = 14)
        vdata <- voomWithQualityWeights(y, design = design, plot = TRUE)
        link_name <- "VoomWithQualityWeightsPlot.pdf"
        link_addr <- "voomwtsplot.pdf"
        link_data <- rbind(link_data, c(link_name, link_addr))
        invisible(dev.off())

        # Generating fit data and top table with weights
        wts <- vdata$weights
        voomfit <- lmFit(vdata, design, weights = wts)
    } else {
        voom_pdf <- make_out("voomplot.pdf")
        voom_png <- make_out("voomplot.png")
        # Creating voom data object and plot
        png(voom_png, width = 500, height = 500)
        vdata <- voom(y, design = design, plot = TRUE)
        img_name <- "VoomPlot"
        img_addr <- "voomplot.png"
        image_data <- rbind(image_data, c(img_name, img_addr))
        invisible(dev.off())

        pdf(voom_pdf)
        vdata <- voom(y, design = design, plot = TRUE)
        link_name <- "VoomPlot.pdf"
        link_addr <- "voomplot.pdf"
        link_data <- rbind(link_data, c(link_name, link_addr))
        invisible(dev.off())

        # Generate voom fit
        voomfit <- lmFit(vdata, design)
    }

    # Save normalised counts (log2cpm)
    if (want_norm) {
        norm_counts <- data.frame(vdata$genes, vdata$E, check.names = FALSE)
        write.table(norm_counts, file = norm_out, row.names = FALSE, sep = "\t", quote = FALSE)
        link_data <- rbind(link_data, c((paste0(de_method, "_", "normcounts.tsv")), (paste0(de_method, "_", "normcounts"))))
    }

    # Fit linear model and estimate dispersion with eBayes
    voomfit <- contrasts.fit(voomfit, contrasts)
    if (want_robust) {
        fit <- eBayes(voomfit, robust = TRUE)
    } else {
        fit <- eBayes(voomfit, robust = FALSE)
    }
    plot_data <- vdata
}

# plot final model mean-variance trend with plotSA
sa_png <- make_out("saplot.png")
sa_pdf <- make_out("saplot.pdf")

png(sa_png, width = 500, height = 500)
plotSA(fit, main = "Final model: Mean-variance trend (SA Plot)")
img_name <- "SAPlot.png"
img_addr <- "saplot.png"
image_data <- rbind(image_data, c(img_name, img_addr))
invisible(dev.off())

pdf(sa_pdf)
plotSA(fit, main = "Final model: Mean-variance trend (SA Plot)")
link_name <- "SAPlot.pdf"
link_addr <- "saplot.pdf"
link_data <- rbind(link_data, c(link_name, link_addr))
invisible(dev.off())

# Save library size info
if (want_libinfo) {
    efflibsize <- round(y$samples$lib.size * y$samples$norm.factors)
    libsizeinfo <- cbind(y$samples, EffectiveLibrarySize = efflibsize)
    libsizeinfo$lib.size <- round(libsizeinfo$lib.size) # nolint
    names(libsizeinfo)[names(libsizeinfo) == "sampleID"] <- "SampleID"
    names(libsizeinfo)[names(libsizeinfo) == "lib.size"] <- "LibrarySize"
    names(libsizeinfo)[names(libsizeinfo) == "norm.factors"] <- "NormalisationFactor"
    write.table(libsizeinfo, file = "libsizeinfo", row.names = FALSE, sep = "\t", quote = FALSE)
}

print("Generating DE results")

if (want_treat) {
    print("Applying TREAT method")
    if (want_robust) {
        fit <- treat(fit, lfc = opt$lfcReq, robust = TRUE)
    } else {
        fit <- treat(fit, lfc = opt$lfcReq, robust = FALSE)
    }
}

status <- decideTests(fit,
    adjust.method = opt$pAdjOpt, p.value = opt$pValReq,
    lfc = opt$lfcReq
)
sum_status <- summary(status)

for (i in seq_along(cons)) {
    con_name <- cons[i]
    con <- cons[i]
    con <- gsub("\\(|\\)", "", con)
    # Collect counts for differential expression
    up_count[i] <- sum_status["Up", i]
    down_count[i] <- sum_status["Down", i]
    flat_count[i] <- sum_status["NotSig", i]

    # Write top expressions table
    if (want_treat) {
        top <- topTreat(fit, coef = i, adjust.method = opt$pAdjOpt, number = Inf, sort.by = "P")
    } else {
        top <- topTable(fit, coef = i, adjust.method = opt$pAdjOpt, number = Inf, sort.by = "P")
    }
    write.table(top, file = top_out[i], row.names = FALSE, sep = "\t", quote = FALSE)
    link_name <- paste0(de_method, "_", con, ".tsv")
    link_addr <- paste0(de_method, "_", con, ".tsv")
    link_data <- rbind(link_data, c(link_name, link_addr))

    # Plot MD (log ratios vs mean average) using limma package on weighted
    pdf(md_pdf[i])
    limma::plotMD(fit,
        status = status[, i], coef = i,
        main = paste("MD Plot:", unmake_names(con)),
        hl.col = alpha(c("firebrick", "blue"), 0.4), values = c(1, -1),
        xlab = "Average Expression", ylab = "logFC"
    )
    abline(h = 0, col = "grey", lty = 2)
    link_name <- paste0("MDPlot_", con, ".pdf")
    link_addr <- paste0("mdplot_", con, ".pdf")
    link_data <- rbind(link_data, c(link_name, link_addr))
    invisible(dev.off())

    # Generate Glimma interactive Volcano, MD plot and tables, requires annotation file (assumes gene labels/symbols in 2nd column)
    if ("i" %in% plots & have_anno) {
        # make gene labels unique to handle NAs
        geneanno <- y$genes
        geneanno[, 2] <- make.unique(geneanno[, 2])

        # use the logcpms for the counts
        if (want_trend) {
            cnts <- logcpm
        } else {
            cnts <- vdata$E
        }

        # MD plot
        Glimma::glMDPlot(fit,
            coef = i, counts = cnts, anno = geneanno, groups = factors[, 1],
            status = status[, i], sample.cols = col.group,
            main = paste("MD Plot:", unmake_names(con)), side.main = colnames(y$genes)[2],
            folder = paste0("glimma_", unmake_names(con)), launch = FALSE
        )
        link_name <- paste0("Glimma_MDPlot_", con, ".html")
        link_addr <- paste0("glimma_", con, "/MD-Plot.html")
        link_data <- rbind(link_data, c(link_name, link_addr))

        # Volcano plot
        Glimma::glXYPlot(
            x = fit$coefficients[, i], y = -log10(fit$p.value[, i]), counts = cnts, anno = geneanno, groups = factors[, 1],
            status = status[, i], sample.cols = col.group,
            main = paste("Volcano Plot:", unmake_names(con)), side.main = colnames(y$genes)[2],
            xlab = "logFC", ylab = "-log10(P-value)",
            folder = paste0("glimma_volcano_", unmake_names(con)), launch = FALSE
        )
        link_name <- paste0("Glimma_VolcanoPlot_", con, ".html")
        link_addr <- paste0("glimma_volcano_", con, "/XY-Plot.html")
        link_data <- rbind(link_data, c(link_name, link_addr))
    }

    # Plot Volcano
    pdf(vol_pdf[i])
    if (have_anno) {
        # labels must be in second column currently
        labels <- fit$genes[, 2]
    } else {
        labels <- fit$genes$GeneID
    }
    limma::volcanoplot(fit,
        coef = i,
        main = paste("Volcano Plot:", unmake_names(con)),
        highlight = opt$topgenes,
        names = labels
    )
    link_name <- paste0("VolcanoPlot_", con, ".pdf")
    link_addr <- paste0("volplot_", con, ".pdf")
    link_data <- rbind(link_data, c(link_name, link_addr))
    invisible(dev.off())

    # PNG of MD and Volcano
    png(mdvol_png[i], width = 1000, height = 500)
    par(mfrow = c(1, 2), mar = c(5, 4, 2, 2) + 0.1, oma = c(0, 0, 3, 0))

    # MD plot
    limma::plotMD(fit,
        status = status[, i], coef = i, main = "MD Plot",
        hl.col = alpha(c("firebrick", "blue"), 0.4), values = c(1, -1),
        xlab = "Average Expression", ylab = "logFC"
    )
    abline(h = 0, col = "grey", lty = 2)

    # Volcano
    if (have_anno) {
        # labels must be in second column currently
        limma::volcanoplot(fit,
            coef = i, main = "Volcano Plot",
            highlight = opt$topgenes,
            names = fit$genes[, 2]
        )
    } else {
        limma::volcanoplot(fit,
            coef = i, main = "Volcano Plot",
            highlight = opt$topgenes,
            names = fit$genes$GeneID
        )
    }

    img_name <- paste0("MDVolPlot_", con)
    img_addr <- paste0("mdvolplot_", con, ".png")
    image_data <- rbind(image_data, c(img_name, img_addr))
    title(paste0("Contrast: ", con_name), outer = TRUE, cex.main = 1.5)
    invisible(dev.off())

    if ("h" %in% plots) {
        # Plot Heatmap
        topgenes <- rownames(top[1:opt$topgenes, ])
        if (want_trend) {
            topexp <- plot_data[topgenes, ]
        } else {
            topexp <- plot_data$E[topgenes, ]
        }
        pdf(heat_pdf[i])
        mycol <- colorpanel(1000, "blue", "white", "red")
        if (have_anno) {
            # labels must be in second column currently
            labels <- top[topgenes, 2]
        } else {
            labels <- rownames(topexp)
        }
        heatmap.2(topexp,
            scale = "row", Colv = FALSE, Rowv = FALSE, dendrogram = "none",
            main = paste("Contrast:", unmake_names(con), "\nTop", opt$topgenes, "genes by adj.P.Val"),
            trace = "none", density.info = "none", lhei = c(2, 10), margin = c(8, 6), labRow = labels, cexRow = 0.7, srtCol = 45,
            col = mycol, ColSideColors = col.group
        )
        link_name <- paste0("Heatmap_", con, ".pdf")
        link_addr <- paste0("heatmap_", con, ".pdf")
        link_data <- rbind(link_data, c(link_name, link_addr))
        invisible(dev.off())
    }

    if ("s" %in% plots) {
        # Plot Stripcharts of top genes
        pdf(strip_pdf[i], title = paste("Contrast:", unmake_names(con)))
        par(mfrow = c(3, 2), cex.main = 0.8, cex.axis = 0.8)
        cols <- unique(col.group)

        for (j in seq_along(topgenes)) {
            lfc <- round(top[topgenes[j], "logFC"], 2)
            pval <- round(top[topgenes[j], "adj.P.Val"], 5)
            if (want_trend) {
                stripchart(plot_data[topgenes[j], ] ~ factors[, 1],
                    vertical = TRUE, las = 2, pch = 16, cex = 0.8, cex.lab = 0.8, col = cols,
                    method = "jitter", ylab = "Normalised log2 expression", main = paste0(labels[j], "\nlogFC=", lfc, ", adj.P.Val=", pval)
                )
            } else {
                stripchart(plot_data$E[topgenes[j], ] ~ factors[, 1],
                    vertical = TRUE, las = 2, pch = 16, cex = 0.8, cex.lab = 0.8, col = cols,
                    method = "jitter", ylab = "Normalised log2 expression", main = paste0(labels[j], "\nlogFC=", lfc, ", adj.P.Val=", pval)
                )
            }
        }
        link_name <- paste0("Stripcharts_", con, ".pdf")
        link_addr <- paste0("stripcharts_", con, ".pdf")
        link_data <- rbind(link_data, c(link_name, link_addr))
        invisible(dev.off())
    }
}
sig_diff <- data.frame(Up = up_count, Flat = flat_count, Down = down_count)
row.names(sig_diff) <- cons

# Save relevant items as rda object
if (want_rda) {
    print("Saving RData")
    if (want_weight) {
        save(counts, data, y, status, plot_data, labels, factors, wts, fit, top, contrast_data, contrasts, design,
            file = rda_out, ascii = TRUE
        )
    } else {
        save(counts, data, y, status, plot_data, labels, factors, fit, top, contrast_data, contrasts, design,
            file = rda_out, ascii = TRUE
        )
    }
    link_data <- rbind(link_data, c((paste0(de_method, "_analysis.RData")), (paste0(de_method, "_analysis.RData"))))
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
cata("<h3>Limma Analysis Output:</h3>\n")
cata("Links to PDF copies of plots are in 'Plots' section below <br />\n")

for (i in seq_len(nrow(image_data))) {
    if (grepl("density|box|mds|mdvol|wts", image_data$Link[i])) {
        html_image(image_data$Link[i], image_data$Label[i], width = 1000)
    } else {
        html_image(image_data$Link[i], image_data$Label[i])
    }
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
# PDFs
for (i in seq_len(nrow(link_data))) {
    if (grepl(".pdf", link_data$Link[i]) & grepl("density|cpm|boxplot|mds|mdplots|voom|saplot", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

for (i in seq_len(nrow(link_data))) {
    if (grepl("mdplot_", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

for (i in seq_len(nrow(link_data))) {
    if (grepl("volplot", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

for (i in seq_len(nrow(link_data))) {
    if (grepl("heatmap", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

for (i in seq_len(nrow(link_data))) {
    if (grepl("stripcharts", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

cata("<h4>Tables:</h4>\n")
for (i in seq_len(nrow(link_data))) {
    if (grepl("counts$", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    } else if (grepl(".tsv", link_data$Link[i])) {
        html_link(link_data$Link[i], link_data$Label[i])
    }
}

if (want_rda) {
    cata("<h4>R Data Object:</h4>\n")
    for (i in seq_len(nrow(link_data))) {
        if (grepl(".RData", link_data$Link[i])) {
            html_link(link_data$Link[i], link_data$Label[i])
        }
    }
}

if ("i" %in% plots) {
    cata("<h4>Glimma Interactive Results:</h4>\n")
    for (i in seq_len(nrow(link_data))) {
        if (grepl("glimma", link_data$Link[i])) {
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
if (want_trend) {
    list_item("The limma-trend method was used.")
} else {
    list_item("The limma-voom method was used.")
}
if (want_weight) {
    list_item("Weights were applied to samples.")
} else {
    list_item("Weights were not applied to samples.")
}
if (want_treat) {
    list_item(paste("Testing significance relative to a fold-change threshold (TREAT) was performed using a threshold of log2 =", opt$lfcReq, "at FDR of", opt$pValReq, "."))
}
if (want_robust) {
    if (want_treat) {
        list_item("TREAT was used with robust settings (robust = TRUE).")
    } else {
        list_item("eBayes was used with robust settings (robust = TRUE).")
    }
}
if (opt$pAdjOpt != "none") {
    if (opt$pAdjOpt == "BH" || opt$pAdjOpt == "BY") {
        temp_str <- paste0(
            "MD Plot highlighted genes are significant at FDR ",
            "of ", opt$pValReq, " and exhibit log2-fold-change of at ",
            "least ", opt$lfcReq, "."
        )
        list_item(temp_str)
    } else if (opt$pAdjOpt == "holm") {
        temp_str <- paste0(
            "MD Plot highlighted genes are significant at adjusted ",
            "p-value of ", opt$pValReq, "  by the Holm(1979) ",
            "method, and exhibit log2-fold-change of at least ",
            opt$lfcReq, "."
        )
        list_item(temp_str)
    }
} else {
    temp_str <- paste0(
        "MD Plot highlighted genes are significant at p-value ",
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

for (i in seq_len(nrow(factors))) {
    cata("<tr>\n")
    table_head_item(row.names(factors)[i])
    for (j in seq_len(ncol(factors))) {
        table_item(as.character(unmake_names(factors[i, j])))
    }
    cata("</tr>\n")
}
cata("</table>")

cit <- character()
link <- character()
link[1] <- paste0(
    "<a href=\"",
    "http://www.bioconductor.org/packages/release/bioc/",
    "vignettes/limma/inst/doc/usersguide.pdf",
    "\">", "limma User's Guide", "</a>."
)

link[2] <- paste0(
    "<a href=\"",
    "http://www.bioconductor.org/packages/release/bioc/",
    "vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf",
    "\">", "edgeR User's Guide", "</a>"
)

cit[1] <- paste("Please cite the following paper for this tool:")

cit[2] <- paste(
    "Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME,",
    "Asselin-Labat ML, Smyth GK, Ritchie ME (2015). Why weight? ",
    "Modelling sample and observational level variability improves power ",
    "in RNA-seq analyses. Nucleic Acids Research, 43(15), e97."
)

cit[3] <- paste(
    "Please cite the paper below for the limma software itself.",
    "Please also try to cite the appropriate methodology articles",
    "that describe the statistical methods implemented in limma,",
    "depending on which limma functions you are using. The",
    "methodology articles are listed in Section 2.1 of the",
    link[1],
    "Cite no. 3 only if sample weights were used."
)
cit[4] <- paste(
    "Smyth GK (2005). Limma: linear models for microarray data.",
    "In: 'Bioinformatics and Computational Biology Solutions using",
    "R and Bioconductor'. R. Gentleman, V. Carey, S. doit,.",
    "Irizarry, W. Huber (eds), Springer, New York, pages 397-420."
)
cit[5] <- paste(
    "Please cite the first paper for the software itself and the",
    "other papers for the various original statistical methods",
    "implemented in edgeR.  See Section 1.2 in the", link[2],
    "for more detail."
)
cit[6] <- paste(
    "Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a",
    "Bioconductor package for differential expression analysis",
    "of digital gene expression data. Bioinformatics 26, 139-140"
)
cit[7] <- paste(
    "Robinson MD and Smyth GK (2007). Moderated statistical tests",
    "for assessing differences in tag abundance. Bioinformatics",
    "23, 2881-2887"
)
cit[8] <- paste(
    "Robinson MD and Smyth GK (2008). Small-sample estimation of",
    "negative binomial dispersion, with applications to SAGE data.",
    "Biostatistics, 9, 321-332"
)
cit[9] <- paste(
    "McCarthy DJ, Chen Y and Smyth GK (2012). Differential",
    "expression analysis of multifactor RNA-Seq experiments with",
    "respect to biological variation. Nucleic Acids Research 40,",
    "4288-4297"
)
cit[10] <- paste(
    "Law CW, Chen Y, Shi W, and Smyth GK (2014). Voom:",
    "precision weights unlock linear model analysis tools for",
    "RNA-seq read counts. Genome Biology 15, R29."
)
cit[11] <- paste(
    "Ritchie ME, Diyagama D, Neilson J, van Laar R,",
    "Dobrovic A, Holloway A and Smyth GK (2006).",
    "Empirical array quality weights for microarray data.",
    "BMC Bioinformatics 7, Article 261."
)
cata("<h3>Citations</h3>\n")
cata(cit[1], "\n")
cata("<br>\n")
cata(cit[2], "\n")

cata("<h4>limma</h4>\n")
cata(cit[3], "\n")
cata("<ol>\n")
list_item(cit[4])
list_item(cit[10])
list_item(cit[11])
cata("</ol>\n")

cata("<h4>edgeR</h4>\n")
cata(cit[5], "\n")
cata("<ol>\n")
list_item(cit[6])
list_item(cit[7])
list_item(cit[8])
list_item(cit[9])
cata("</ol>\n")

cata("<p>Please report problems or suggestions to: su.s@wehi.edu.au</p>\n")

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
