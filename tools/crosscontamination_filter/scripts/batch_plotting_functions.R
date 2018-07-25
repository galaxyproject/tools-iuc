#!/usr/bin/env R

##
## Batch Plotting Functions
##
calculatePlateIndexes <- function(plate.form, full.barcode.size, num.plates){
    #' Determine plotting positions of plate lines (under false model)
    #'
    #' Assumes all plates are the same size and span the full range of the
    #' barcodes.
    #'
    #' @param plate.form list of vectors mapping plates to batches
    #' @param full.barcode.size size of the complete barcodes list
    #' @param num.plates, number of plates
    #' @return sequence of discrete plate-boundary positions
    batches.per.plate <- length(plate.form[[1]])
    plate.size <- batches.per.plate * full.barcode.size

    return(seq(0, num.plates * plate.size, plate.size))
}


calculateFullBarcodeIndexes <- function(num.batches, full.barcode.size){
    #' Determines plotting position of batch lines (under false model)
    #'
    #' For N batches and a list of actually detected barcodes in the header,
    #' generates where the blue lines should be
    #'
    #' @param num.batches number of batches in experiment
    #' @param full.barcode.size size of all barcodes
    #' @return sequence of discrete batch positions
    bsize <- full.barcode.size
    return(seq(0, num.batches * bsize, bsize))
}



calculateRealBarcodeIndexes <- function(barcode.form, full.barcode.size){
    #' Determine plotting position of the true batch lines (under true model)
    #'
    #' For N batches a list of actually USED barcodes as given by the spec,
    #' generates where the green lines should be
    #'
    #' @param barcode.form list of barcode formats and the batches they map to
    #' @param full.barcode.size size of all barcodes
    #' @return list of useful vectors: true batch positions using whole matrix,
    #'         true batch positions using the filtered matrix which contains
    #'         only real barcodes, and a list of batches and their respective
    #'         sizes.
    batches <- c()
    res <- sapply(names(barcode.form), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        size.of.range <- length(seq(rng[1],rng[2]))
        sub.batches <- barcode.form[[key]]  # 1,3,5,7 or 2,4,6,8
        res2 <- lapply(sub.batches, function(bat){
            batches[[bat]] <<- size.of.range
        })
    })
    ## We now have sizes per batch, in order of batch
    ## Need to place these at positions after each full barcode size
    positions <- c()
    real_positions <- c(0)

    res <- sapply(1:length(batches), function(b){
        batch.start <- (b-1) * full.barcode.size
        batch.size <- batches[[b]]
        positions <<- c(positions, batch.start + batch.size)
        real_positions[[b+1]] <<- batch.size + real_positions[[b]]
    })

    real_positions <- real_positions[2:length(real_positions)]   
    
    return(list(unfiltered=positions,filtered=real_positions, batches=batches))
}

calculateRealPlateIndexes <- function(plate.form, batches, num.plates){
    #' Determine true plate positions given variable batch sizes
    #'
    #' Given the true size of each batch, and which batches exist in which plates
    #' calculate the size of each plate
    #'
    #' @param plate.form list of vectors mapping plates to batches
    #' @param batches list of batches and their respective sizes
    #' @param num.plates number of plates
    #' @return sequence of plate positions
    batches.per.plate <- length(plate.form[[1]])

    size.of.plate <- 0
    res <- sapply(plate.form[[1]], function(batch){
        batch.size <- batches[[batch]]
        size.of.plate <<- size.of.plate + batch.size
    })

    return(seq(0, num.plates * size.of.plate, size.of.plate))
}
