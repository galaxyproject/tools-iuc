## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("DEXSeq")
    library('getopt')
})

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'rdata', 'r', 1, "character",
    'primaryfactor', 'p', 1, "character",
    'geneid', 'g', 1, "character",
    'fdr', 'c', 1, "double",
    'transcripts', 't', 1, "logical",
    'names', 'a', 1, "logical",
    'normcounts', 'n', 1, "logical",
    'splicing', 's', 1, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);  

res <- readRDS(opt$rdata)

pdf("plot.pdf")
plotDEXSeq(res, opt$geneid, FDR=opt$fdr, fitExpToVar=opt$primaryfactor,
    norCounts=opt$normcounts, expression=TRUE, splicing=opt$splicing,
    displayTranscripts=opt$transcripts, names=opt$names, legend=TRUE,
    color=NULL, color.samples=NULL, transcriptDb=NULL)
dev.off()

sessionInfo()