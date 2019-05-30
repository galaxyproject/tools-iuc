#!/usr/bin/env R
VERSION = "0.2"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

source(args[1])

## debug in debug dir
source('../debug/config.debug')

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


barcode.data <- calculateBarcodePositions(spec$format, num.barcodes)
plate.data <- calculatePlatePositions(spec$plates, num.barcodes, barcode.data)


ordering <- reorderMatrixHeaders(barcodes, colnames(input_matrix), spec$format, spec$plates)

## Unfiltered, but sorted matrix
nmatrix <- input_matrix[,ordering$all]
plot.prefilter <- contaminationPlot("Pre-Filter", colSums(nmatrix),
                                    barcode.data, plate.data, RAW=TRUE)

## Filtered, but sorted matrix
cmatrix <- input_matrix[,ordering$filtered]
plot.postfilter <- contaminationPlot("Post-Filter", colSums(cmatrix),
                                     barcode.data, plate.data, RAW=FALSE)

plot.histogram.pre <- log10histoPlot("Histogram of Pre-Filter Matrix Counts", colSums(nmatrix))
plot.histogram.post <- log10histoPlot("Histogram of Post-Filter Matrix Counts", colSums(cmatrix))

pdf(out.pdf)
plot.prefilter
plot.postfilter
plot.histogram.pre
plot.histogram.post
dev.off()

write.table(cmatrix, file=out.table, quote=FALSE, na="0", sep="\t")
