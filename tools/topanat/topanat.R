# This tool performs a GO-like enrichment of anatomical terms and cell types
# mapped to genes by expression patterns, using the BgeeDB R package (topAnat).
#
# ARGS: species",        "s", 1, "character"  -Species name in the format Genus_species (e.g. Homo_sapiens)
#       dataTypes",      "d", 2, "character"  -Comma-separated datatypes: rna_seq, sc_full_length, sc_droplet_based, affymetrix
#       stageId",        "t", 2, "character"  -Developmental stage ID to filter expression data (e.g. UBERON:0000092)
#       fgGenes",        "g", 2, "character"  -Newline-separated list of foreground gene IDs (mutually exclusive with fgFile)
#       fgFile",         "f", 2, "character"  -Path to one-column TSV file (no header) containing foreground gene IDs (mutually exclusive with fgGenes)
#       bgGenes",        "G", 2, "character"  -Newline-separated list of background gene IDs (mutually exclusive with bgFile)
#       bgFile",         "b", 2, "character"  -Path to one-column TSV file (no header) containing background gene IDs (mutually exclusive with bgGenes)
#       algorithm",      "a", 2, "character"  -Decorrelation algorithm (default: classic). See topGO::whichAlgorithms()
#       statistics",     "S", 2, "character"  -Test statistic (default: fisher). See topGO::whichTests()
#       nodeSize",       "n", 2, "integer"    -Minimum number of genes mapped to a node for it to be tested (default: 10)
#       confidence",     "c", 2, "character"  -Call quality: silver (default), gold, all, high_quality
#       resultFile",     "r", 1, "character"  -Path to output TSV file with enrichment results
#
# OUT:
#       TSV table with GO-like enrichment results for anatomical entities / cell types
#
# Author: based on edgeR Galaxy wrapper structure

# Record starting time
time_start <- as.character(Sys.time())

# Setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
})

# We need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Load all required libraries
library(methods, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)
library(BgeeDB, quietly = TRUE, warn.conflicts = FALSE)

################################################################################
### Input Processing
################################################################################

# Collect arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec defined by the enclosed list.
spec <- matrix(
    c(
        "species",    "s", 1, "character",
        "dataTypes",  "d", 2, "character",
        "stageId",    "t", 2, "character",
        "fgGenes",    "g", 2, "character",
        "fgFile",     "f", 2, "character",
        "bgGenes",    "G", 2, "character",
        "bgFile",     "b", 2, "character",
        "algorithm",  "a", 2, "character",
        "statistics", "S", 2, "character",
        "nodeSize",   "n", 2, "integer",
        "confidence", "c", 2, "character",
        "resultFile", "r", 1, "character"
    ),
    byrow = TRUE, ncol = 4
)
opt <- getopt(spec)

# Validate required arguments
if (is.null(opt$species)) {
    cat("A species name is required (e.g. Homo_sapiens).\n", file = stderr())
    q(status = 1)
}

if (is.null(opt$fgGenes) && is.null(opt$fgFile)) {
    cat("Either foreground genes (--fgGenes) or a foreground genes file (--fgFile) is required.\n", file = stderr())
    q(status = 1)
}

if (is.null(opt$resultFile)) {
    cat("An output result file path is required.\n", file = stderr())
    q(status = 1)
}

# Set defaults for optional arguments
if (is.null(opt$algorithm)) {
    opt$algorithm <- "classic"
}

if (is.null(opt$statistics)) {
    opt$statistics <- "fisher"
}

if (is.null(opt$nodeSize)) {
    opt$nodeSize <- 10L
}

if (is.null(opt$confidence)) {
    opt$confidence <- "silver"
}

# Parse dataTypes: accept a comma-separated string and split into a character vector
if (is.null(opt$dataTypes) || opt$dataTypes == "") {
    dataTypes <- character(0)
} else {
    dataTypes <- unlist(strsplit(opt$dataTypes, split = ",", fixed = TRUE))
    dataTypes <- trimws(dataTypes)
}

# Parse foreground genes: accept a newline-separated string (from a Galaxy text area) and split into a character vector
if (!is.null(opt$fgGenes)) {
    foregroundGenes <- unlist(strsplit(opt$fgGenes, split = "\n", fixed = TRUE))
    foregroundGenes <- trimws(foregroundGenes)
    foregroundGenes <- foregroundGenes[nzchar(foregroundGenes)]
} else {
    foregroundGenes <- NULL
}

# Parse background genes: accept a newline-separated string (from a Galaxy text area) and split into a character vector
if (!is.null(opt$bgGenes)) {
    backgroundGenes <- unlist(strsplit(opt$bgGenes, split = "\n", fixed = TRUE))
    backgroundGenes <- trimws(backgroundGenes)
    backgroundGenes <- backgroundGenes[nzchar(backgroundGenes)]
} else {
    backgroundGenes <- NULL
}

################################################################################
### Run topAnat
################################################################################

topAnat_galaxy(
    species             = opt$species,
    dataTypes           = dataTypes,
    stageId             = opt$stageId,
    foregroundGenes     = foregroundGenes,
    backgroundGenes     = backgroundGenes,
    foregroundGenesFile = opt$fgFile,
    backgroundGenesFile = opt$bgFile,
    algorithm           = opt$algorithm,
    statistics          = opt$statistics,
    nodeSize            = opt$nodeSize,
    confidence          = opt$confidence,
    resultFile          = opt$resultFile
)

cat("TopAnat analysis completed successfully.\n")
cat(paste("Started:", time_start, "\n"))
cat(paste("Finished:", as.character(Sys.time()), "\n"))

time_end <- as.character(Sys.time())
time_taken <- capture.output(round(difftime(time_end, time_start), digits = 3))
time_taken <- gsub("Time difference of ", "", time_taken, fixed = TRUE)
cat(paste("Time taken:", time_taken, "seconds\n"))
