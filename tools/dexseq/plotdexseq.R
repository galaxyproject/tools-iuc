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
})

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    "rdata", "r", 1, "character",
    "primaryfactor", "p", 1, "character",
    "geneid", "g", 1, "character",
    "genefile", "f", 1, "character",
    "fdr", "c", 1, "double",
    "transcripts", "t", 1, "logical",
    "names", "a", 1, "logical",
    "normcounts", "n", 1, "logical",
    "splicing", "s", 1, "logical",
    "pl_width", "w", 2, "integer",
    "pl_height", "h", 2, "integer"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

res <- readRDS(opt$rdata)

if (!is.null(opt$genefile)) {
    genes <- read.delim(opt$genefile, header = FALSE)
    genes <- genes[, 1]
} else {
    genes <- opt$geneid
}

pl_width <- pl_height <- 7
if (!is.null(opt$pl_width)) pl_width <- opt$pl_width
if (!is.null(opt$pl_height)) pl_height <- opt$pl_height
pdf("plot.pdf", width = pl_width, height = pl_height)
for (i in genes) {
    par(oma = c(pl_height * 0.2, pl_width * 0.2, pl_height * 0.2, pl_width * 0.2))
    plotDEXSeq(res, i,
        FDR = opt$fdr, fitExpToVar = opt$primaryfactor,
        norCounts = opt$normcounts, expression = TRUE, splicing = opt$splicing,
        displayTranscripts = opt$transcripts, names = opt$names, legend = TRUE,
        color = NULL, color.samples = NULL, transcriptDb = NULL
    )
}
dev.off()

sessionInfo()
