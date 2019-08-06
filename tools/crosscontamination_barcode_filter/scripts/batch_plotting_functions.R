#!/usr/bin/env R

suppressPackageStartupMessages(require(data.table))

##
## Batch Plotting Functions
##
sortedBatchesOrPlates <- function(batch.list){
    #' B0, B1, B11, B12, B2, B3, ...
    #' to B0,B1,B2,---B11,B12
    vals.and.index = sort(as.integer(sub("^[BP]", "", batch.list)), index.return=TRUE)
    return(batch.list[vals.and.index$ix])
}


calculateBarcodePositions <- function(barcode.form, full.barcode.size){
    #' Determine x-axis positions of all batches under the context of
    #' unfiltered barcodes (full set), filtered (real set), dividing line
    #' (the position of real set in the full set).
    #'
    #' @param barcode.form list of barcode formats and the batches they map to
    #' @param full.barcode.size size of all barcodes
    #' @return dataframe of batch information: sizes, unfiltered, filtered,
    #'         and dividing line.
    sizes <- list(B0=0)

    res <- sapply(names(barcode.form), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        size.of.range <- length(seq(rng[1],rng[2]))
        sub.batches <- barcode.form[[key]]  # 1,3,5,7 or 2,4,6,8
        res2 <- lapply(sub.batches, function(bat){
            sizes[[paste("B",bat, sep="")]] <<- size.of.range
        })
    })

    ## We now have sizes per batch, in order of batch
    ## Need to place these at positions after each full barcode size

    ## Below we have "positions" which has the END positions of each batch under
    ## the assumption of using full barcodes. The "real_positions" contains the
    ## END positions of each batch under the assumption of using only the real
    ## subsetted barcodes.
    unfilter_positions <- list(B0=0)
    filtered_positions <- list(B0=0)
    filter_in_unfilter <- list(B0=0) ## dividing line between real and false barcodes in each batch

    res <- sapply(sortedBatchesOrPlates(names(sizes)), function(batch.name){

        batch.num <- as.integer(sub("B","", batch.name))
        if (batch.num > 0){
            batch.size <- sizes[[batch.name]]  ## 96
            batch.name.previous = paste("B", batch.num-1, sep="")
            batch.start <- unfilter_positions[[batch.name.previous]]
            filt.batch.start <- filtered_positions[[batch.name.previous]]

            unfilter_positions[[batch.name]] <<- batch.start + full.barcode.size
            filtered_positions[[batch.name]] <<- filt.batch.start + batch.size
            filter_in_unfilter[[batch.name]] <<- batch.start + batch.size
        }
    })

    # Put into a dataframe, merging lists on their common names
    dd <- data.frame(rbindlist(list(
        unfilter_positions=unfilter_positions,
        filter_in_unfilter=filter_in_unfilter,
        filtered_positions=filtered_positions,
        sizes=sizes),  ## sizes go last to not mess up the column name ordering
        use.names = TRUE, idcol = TRUE))

    rownames(dd) <- dd$.id
    dd <- dd[,!(colnames(dd) %in% ".id")]

    return(dd)
}

calculatePlatePositions <- function(plate.form, full.barcode.size, all.batch.data){
    #' Determine the x-axis plate positions for each of the unfiltered and filtered sets
    #'
    #' Given the true size of each batch, and which batches exist in which plates
    #' calculate the size of each plate
    #'
    #' @param plate.form list of vectors mapping plates to batches
    #' @param full.barcode.size size of the full set of barcodes
    #' @param all.batch.data the output of 'calculateBarcodePositions'
    #' @return dataframe of plate information pertaining to positions and sizes of plates
    unfilter.plates = list(P0=0)
    filtered.plates = list(P0=0)
    unfilter.plates.sizes = list(P0=0)
    filtered.plates.sizes = list(P0=0)

    res <- sapply(sortedBatchesOrPlates(names(plate.form)), function(plate.num){

        unfilter.plate.size = 0
        filtered.plate.size = 0

        batches <- plate.form[[plate.num]]

        res2 <- sapply(sort(batches), function(batch.num){
            batch.size <- all.batch.data["sizes",paste("B", batch.num, sep="")]
            unfilter.plate.size <<- unfilter.plate.size + full.barcode.size
            filtered.plate.size <<- filtered.plate.size + batch.size
        })

        plate.name = paste("P", plate.num, sep="")
        plate.name.previous = paste("P", as.integer(plate.num) - 1, sep="")

        unfilter.plates.sizes[[plate.name]] <<- unfilter.plate.size
        filtered.plates.sizes[[plate.name]] <<- filtered.plate.size

        filtered.plates[[plate.name]] <<- filtered.plates[[plate.name.previous]] + filtered.plate.size
        unfilter.plates[[plate.name]] <<- unfilter.plates[[plate.name.previous]] + unfilter.plate.size
    })

    # Put into a dataframe, merging lists on their common names
    dd <- data.frame(rbindlist(list(
        unfilter.plates=unfilter.plates,
        unfilter.plates.sizes=unfilter.plates.sizes,
        filtered.plates=filtered.plates,
        filtered.plates.sizes=filtered.plates.sizes),
        use.names = TRUE, idcol = TRUE))

    rownames(dd) <- dd$.id
    dd <- dd[,!(colnames(dd) %in% ".id")]

    return(dd)
}
