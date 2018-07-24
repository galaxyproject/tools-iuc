#!/usr/bin/env R
VERSION = "0.1"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION: ", VERSION))
     stop("Please provide the config file")
}


##                                         # Config file
## input_matrix <- read.table(
##     '../test-data/input_matrix.tsv',
##     stringsAsFactors = F,
##     na.strings=c("NA", "-", "?", "."),
##     header=TRUE,
##     row.names=1
## )
## input_matrix[is.na(input_matrix)] <- 0

## spec = list(
##     barcodes = '../test-data/celseq_barcodes.192.raw',
##     format = list(
##         "85-96"   = c(1,3,5,7),
##         "150-192" = c(2,4,6,8)
##     ),
##     plates = list(
##         "1" = c(1,2,3,4),
##         "2" = c(5,6,7,8)
##     )
## )
## regex.extract = ".*P(\\d)_(\\d)_([ACTG]+)"
## regex.display = "P\\1_B\\2_\\3"
## script.dir = "../scripts"

##                                         # End config file

source(args[1])


convertHeadersToSensible <- function(regex.from, regex.to, matrix){
    return(sub(regex.from, regex.to, colnames(matrix)))
}

colnames(input_matrix) <- convertHeadersToSensible(regex.extract, regex.display, input_matrix)


## At this point we expect our headers to have names: P1_B2_AACCTT

barcodes <- NULL
num.barcodes <- NULL
num.batches <- NULL
num.plates <- NULL

source(file.path(script.dir, "config_assertions.R"))

sanityCheck <- function(spec, matrix.headers){
    barcodes <<- scan(spec$barcodes, what="", sep="\n")
    num.barcodes <<- length(barcodes)

    used.barcodes <- checkNoMissingRanges(spec$format, barcodes)
    assertNoMissingBarcodes(matrix.headers, used.barcodes)
    num.batches <<- assertNoMissingBatches(spec$format, spec$plates)
}

sanityCheck(spec, colnames(input_matrix))

source(file.path(script.dir, "batch_plotting_functions.R"))
source(file.path(script.dir, "reorder_matrix_headers.R"))
source(file.path(script.dir, "contamination_plot.R"))

ordering <- reorderMatrixHeaders(barcodes, colnames(input_matrix), spec$format)


real.indexes = calculateRealBarcodeIndexes(spec$format, num.barcodes)
plate.indexes = calculatePlateIndexes(spec$plates, num.barcodes)
plate.indexes.real = calculateRealPlateIndexes(spec$plates, real.indexes$batches)


## Unfiltered
nmatrix <- input_matrix[,ordering$all]

plot.prefilter <- contaminationPlot(
    colSums(nmatrix),
    title="Pre-Filter",
    plate.indexes,
    calculateFullBarcodeIndexes(num.batches, num.barcodes),
    real.indexes$unfiltered
)


## Filtered
cmatrix <- input_matrix[,ordering$correct]

plot.postfilter <- contaminationPlot(
    colSums(cmatrix),
    title="Post-Filter",
    plate.indexes.real,
    calculateFullBarcodeIndexes(num.batches, num.barcodes),
    real.indexes$filtered,
    filtered = T
)

plot.histogram <- log10histoPlot(colSums(cmatrix))

pdf(out.pdf)
plot.prefilter
plot.postfilter
plot.histogram
dev.off()


write.table(cmatrix, file=out.table, quote=FALSE, na="0", sep="\t")
