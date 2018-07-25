#!/usr/bin/env R
##
## Sanity Check Methods
##
checkNoMissingRanges <- function(format, barcodes){
    #' Checks coverage of barcodes across all specified ranges
    #'
    #' e.g. 1-50, 60-80 -- barcodes 51-59 are not specified. This
    #'       is not a fatal error, but a warning is issued.
    #'
    #' @param format barcode range and the batches they map to
    #' @param barcodes full list of barcodes
    #' @return list of specified barcodes
    ranges <- c()
    res <- sapply(names(format), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        ranges <<- c(ranges, seq(rng[1],rng[2]))
    })
    full.range <- seq(min(ranges),max(ranges))
    not.in <- !(full.range %in% ranges)
    if (sum(not.in) != 0) {
        message("Warning: values[",
                paste(full.range[not.in], collapse=","),
                "] -> barcodes[",
                paste(barcodes[not.in], collapse=","),
                "] are not used."
                )
    } else {
        message("All barcodes accounted for.")
    }
    return(barcodes[!not.in])
}

checkNoMissingBarcodes <- function(headers, barcodes){
    #' Extracts barcodes in the headers and compares them with those in barcodes
    #'
    #' @param headers matrix headers, must be of P1_B2_ACTG format
    #' @param barcodes full list of barcodes
    barcs.in.matrix <- unique(sort(sub(".*_.*_([ACTG]+)", "\\1", headers)))
    not.in <- !(barcs.in.matrix %in% barcodes)
    if (sum(not.in) > 0){
        message("Warning: Barcodes in matrix not in barcodes file\n", barc.in.matrix[not.in])
    } else {
        message("All input matrix barcodes accounted for.")
    }
}

assertNoMissingBatches <- function(format, plates){
    #' Checks the barcode and plate spec match
    #' 
    #' These must specify the same batches.
    #'
    #' @param format barcode format, ranges to batches
    #' @param plates plate format, plates to batches
    #' @return number of batches
    batches.form = c()
    batches.plate = c()
    for (form in format){batches.form = c(batches.form, form)}
    for (plate in plates){batches.plate = c(batches.plate, plate)}

    if (length(batches.plate) != length(batches.form)){
        stop("Error: The number of batches specified in the plate do not match those given in the barcode format")
    }

    range.form <- seq(min(batches.form), max(batches.form))
    range.plate <- seq(min(batches.plate), max(batches.plate))

    if (sum(!(range.form %in% batches.form)) > 0){
        stop("Error: Missing batch in barcode format")
    }
    if (sum(!(range.plate %in% batches.plate)) > 0){
        stop("Error: Missing batch in plate format")
    }
    return(length(range.form))
}

sanityCheck <- function(spec, matrix.headers){
    #' Checks specification and matrix headers for consistency
    #'
    #' @param spec experiment specification
    #' @param matrix.headers column names of input matrix
    #' @return list of barcodes, as well as number of barcodes, plates, and batches
    barcodes <- scan(spec$barcodes, what="", sep="\n")
    num.barcodes <- length(barcodes)
    num.plates <- length(names(spec$plate))
    used.barcodes <- checkNoMissingRanges(spec$format, barcodes)
    num.batches <- assertNoMissingBatches(spec$format, spec$plates)
    checkNoMissingBarcodes(matrix.headers, used.barcodes)

    return(list(barc=barcodes, barc.n=num.barcodes, plates.n=num.plates, batch.n=num.batches))
}
