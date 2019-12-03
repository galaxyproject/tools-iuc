# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "help", "h", 0, "logical",
  "out_file", "o", 1, "character",
  "countsFiles", "n", 1, "character",
  "countsFromAbundance", "r", 1, "character",
  "format", "v", 1, "character",
  "gff_file", "H", 0, "character",
  "tx2gene", "f", 0, "character",
  "geneIdCol", "l", 0, "character",
  "txIdCol" , "p", 1, "character",
  "abundanceCol", "i", 0, "character",
  "countsCol", "y", 1, "character",
  "lengthCol", "x", 1, "character"),
  byrow=TRUE, ncol=4)

opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
	
if (is.null(opt$gff_file) & is.null(opt$tx2gene)) {
  cat("A GFF/GTF file or a tx2gene table is required\n")
  q(status=1)
}

if (opt$format == 'none'){  #custom format
    if (is.null(opt$txIdCol) | is.null(opt$abundanceCol) | is.null(opt$countsCol) | is.null(opt$lengthCol)) {
        cat("If you select a custom format for the input files you need to specify the column names\n")
        q(status=1)
   }
}

if (is.null(opt$countsFiles)) {
  cat("'countsFiles' is required\n")
  q(status=1)
}


# load samples from tab file
samples_df <- read.table(opt$countsFiles, sep="\t", header=TRUE)
colnames(samples_df) <- c("id","path")
rownames(samples_df) <- NULL
# Prepare char vector with files and sample names 
files <- file.path(samples_df[,"path"])
names(files) <- samples_df[,"id"]



library(tximport)

### if the input is a gff/gtf file first need to create the tx2gene table
if (!is.null(opt$gff_file)) {
    suppressPackageStartupMessages({
        library("GenomicFeatures")
    })
    txdb <- makeTxDbFromGFF(opt$gff_file)
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, keys=k, columns="GENEID", keytype="TXNAME")
    # Remove 'transcript:' from transcript IDs (when gffFile is a GFF3 from Ensembl and the transcript does not have a Name)
    tx2gene$TXNAME <- sub('^transcript:', '', tx2gene$TXNAME)

} else {
        tx2gene <- read.table(opt$tx2gene,header=FALSE)
    }



##
if (is.null(opt$geneIdCol)) { ## there is a tx2gene table
    if (opt$format == 'none'){  #predefined format 
        txi_out <- tximport(files, type="none",txIdCol=opt$txIdCol,abundanceCol=opt$abundanceCol,countsCol=opt$countsCol,lengthCol=opt$lengthCol,tx2gene=tx2gene,countsFromAbundance=opt$countsFromAbundance)
    } else {
        txi_out <- tximport(files, type=opt$format, tx2gene=tx2gene,countsFromAbundance=opt$countsFromAbundance)
    }
} else {  # the gene_ID is a column in the counts table
    if (opt$format == 'none'){  #predefined format
        txi_out <- tximport(files, type="none",geneIdCol=opt$geneIdCol,txIdCol=opt$txIdCol,abundanceCol=opt$abundanceCol,countsCol=opt$countsCol,lengthCol=opt$lengthCol,tx2gene=tx2gene,countsFromAbundance=opt$countsFromAbundance)
    } else {
        txi_out <- tximport(files, type=opt$format, geneIdCol=opt$geneIdCol,countsFromAbundance=opt$countsFromAbundance)
    }

}
# write count as table
write.table(txi_out$counts, file=opt$out_file, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
