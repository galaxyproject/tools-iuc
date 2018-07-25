#!/usr/bin/env R
VERSION = "0.1"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION: ", VERSION))
     stop("Please provide the config file")
}

source(args[1])
source(file.path(script.dir, "config_assertions.R"))
source(file.path(script.dir, "batch_plotting_functions.R"))
source(file.path(script.dir, "reorder_matrix_headers.R"))
source(file.path(script.dir, "contamination_plot.R"))


colnames(input_matrix) <- convertHeadersToSensible(
    regex.extract,
    regex.display,
    colnames(input_matrix)
)

sc <- sanityCheck(spec, colnames(input_matrix))

barcodes <- sc$barc
num.barcodes <- sc$barc.n
num.batches <- sc$batch.n
num.plates <- sc$plates.n

real.indexes = calculateRealBarcodeIndexes(spec$format, num.barcodes)
plate.indexes = calculatePlateIndexes(spec$plates, num.barcodes, num.plates)
plate.indexes.real = calculateRealPlateIndexes(spec$plates, real.indexes$batches, num.plates)

ordering <- reorderMatrixHeaders(barcodes, colnames(input_matrix), spec$format)

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

plot.histogram <- log10histoPlot(
    colSums(cmatrix),
    "Histogram of Post-Filter Matrix Counts"
)

pdf(out.pdf)
plot.prefilter
plot.postfilter
plot.histogram
dev.off()


write.table(cmatrix, file=out.table, quote=FALSE, na="0", sep="\t")
