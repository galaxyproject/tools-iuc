#!/usr/bin/env R
                                        #
                                        #
                                        # Batch Plotting Functions
                                        #
calculatePlateIndexes <- function(plate.form, full.barcode.size){
    num.plates <<- length(names(plate.form))
    batches.per.plate <- length(plate.form[[1]])
    plate.size <- batches.per.plate * full.barcode.size

    return(seq(0, num.plates * plate.size, plate.size))
}


calculateFullBarcodeIndexes <- function(num.batches, full.barcode.size){
    #' For N batches and a list of actually detected barcodes in the header,
    #' generates where the blue lines should be
    bsize <- full.barcode.size
    return(seq(0, num.batches * bsize, bsize))
}



calculateRealBarcodeIndexes <- function(barcode.form, full.barcode.size){
    #' For N batches a list of actually USED barcodes as given by the spec,
    #' generates where the green lines should be
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

calculateRealPlateIndexes <- function(plate.form, batches){
    #' Given the true size of each batch, and which batches exist in which plates
    #' calculate the size of each plate
    num.plates <<- length(names(plate.form))
    batches.per.plate <- length(plate.form[[1]])

    size.of.plate <- 0
    res <- sapply(plate.form[[1]], function(batch){
        batch.size <- batches[[batch]]
        size.of.plate <<- size.of.plate + batch.size
    })

    return(seq(0, num.plates * size.of.plate, size.of.plate))
}
