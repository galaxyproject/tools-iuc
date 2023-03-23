## Setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
})
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("DEXSeq")
    library("getopt")
    library("rjson")
})


options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    "verbose", "v", 2, "integer",
    "help", "h", 0, "logical",
    "gtf", "a", 1, "character",
    "outfile", "o", 1, "character",
    "reportdir", "r", 1, "character",
    "rds", "d", 1, "character",
    "factors", "f", 1, "character",
    "threads", "p", 1, "integer",
    "fdr", "c", 1, "double"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

trim <- function(x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

parser <- newJSONParser()
parser$addData(opt$factors)
factors_list <- parser$getObject()

sample_table <- data.frame()
count_files <- c()
factor_names <- c()
primary_factor <- ""
for (factor in factors_list) {
    factor_name <- factor[[1]]
    factor_names <- append(factor_names, paste(factor_name, "exon", sep = ":"))
    factor_values_map_list <- factor[[2]]
    c <- length(factor_values_map_list)
    for (factorValuesMap in factor_values_map_list) {
        for (files in factorValuesMap) {
            for (file in files) {
                if (primary_factor == "") {
                    count_files <- append(count_files, file)
                }
                sample_table[basename(file), factor_name] <- paste(c, names(factorValuesMap), sep = "_")
            }
        }
        c <- c - 1
    }
    if (primary_factor == "") {
        primary_factor <- factor_name
    }
}

factor_names <- append(factor_names, "exon")
factor_names <- append(factor_names, "sample")
factor_names <- rev(factor_names)
formula_full_model <- as.formula(paste("", paste(factor_names, collapse = " + "), sep = " ~ "))
factor_names <- head(factor_names, -1)
formula_reduced_model <- as.formula(paste("", paste(factor_names, collapse = " + "), sep = " ~ "))

sample_table
formula_full_model
formula_reduced_model
primary_factor
count_files
opt$reportdir
opt$threads
getwd()

dxd <- DEXSeqDataSetFromHTSeq(count_files, sampleData = sample_table, design = formula_full_model, flattenedfile = opt$gtf)

colData(dxd)
dxd <- estimateSizeFactors(dxd)
print("Estimated size factors")
sizeFactors(dxd)
bpparam <- MulticoreParam(workers = opt$threads)
dxd <- estimateDispersions(dxd, formula = formula_full_model, BPPARAM = bpparam)
print("Estimated dispersions")
dxd <- testForDEU(dxd, reducedModel = formula_reduced_model, fullModel = formula_full_model, BPPARAM = bpparam)
print("tested for DEU")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = primary_factor, BPPARAM = bpparam)
print("Estimated fold changes")
res <- DEXSeqResults(dxd)
print("Results")
table(res$padj <= opt$fdr)
res_sorted <- res[order(res$padj), ]
head(res_sorted)

export_table <- as.data.frame(res_sorted)
last_column <- ncol(export_table)
for (i in seq_len(nrow(export_table))) {
  export_table[i, last_column] <- paste(export_table[i, last_column][[1]], collapse = ", ")
}
export_table[, c(last_column)] <- sapply(export_table[, c(last_column)], as.character)
write.table(export_table, file = opt$outfile, sep = "\t", quote = FALSE, col.names = FALSE)
print("Written Results")

if (!is.null(opt$rds)) {
    saveRDS(res, file = "DEXSeqResults.rds")
}

if (!is.null(opt$reportdir)) {
    DEXSeqHTML(res, fitExpToVar = primary_factor, path = opt$reportdir, FDR = opt$fdr, color = c("#B7FEA0", "#FF8F43", "#637EE9", "#FF0000", "#F1E7A1", "#C3EEE7", "#CEAEFF", "#EDC3C5", "#AAA8AA"))
}
sessionInfo()
