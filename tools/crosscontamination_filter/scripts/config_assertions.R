#!/usr/bin/env R
                                        #
                                        #
                                        # Sanity Check Methods
                                        #
checkNoMissingRanges <- function(format, barcodes){
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

assertNoMissingBarcodes <- function(headers, barcodes){
    ## headers are in format P1_B1_ACTG
    barcs.in.matrix <- unique(sort(sub(".*_.*_([ACTG]+)", "\\1", headers)))
    not.in <- !(barcs.in.matrix %in% barcodes)
    if (sum(not.in) > 0){
        message("Error: Barcodes in matrix not in barcodes file\n", barc.in.matrix[not.in])
        stop("")
    } else {
        message("All input matrix barcodes accounted for.")
    }
}

assertNoMissingBatches <- function(format, plates){
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
