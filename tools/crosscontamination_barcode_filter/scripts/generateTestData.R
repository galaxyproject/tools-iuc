num.genes = 500
num.cells.per.batch.good = 96
num.cells.per.batch.bad = 11
num.cells.per.batch.empty = 85

barcodes = read.table('../test-data/celseq_barcodes.192.raw', stringsAsFactors=F)[,1]

barcodes.1_96 <- barcodes[1:96]
barcodes.97_192 <- barcodes[97:192]


goodCell <- function(){
    return(rnbinom(n=num.genes, size=0.1, mu=6))
}

badCell <- function(){
    return(rnbinom(n=num.genes, size=0.1, mu=0.5))
}

emptyCell <- function(){
    return(rnbinom(n=num.genes, size=0.1, mu=0.1))
}

generateGoodBatch <- function(batch.prefix="P1_B1", bar.codes){
    tmp.cell <- goodCell()
    for (c in 2:num.cells.per.batch.good){
        tmp.cell <- cbind(tmp.cell, goodCell())
    }
    bar.good <- paste(batch.prefix, bar.codes, sep="_")
    colnames(tmp.cell) <- bar.good
    return(tmp.cell)
}


generateBadBatch <- function(batch.prefix="P1_B1", bar.codes){
    tmp.cell <- emptyCell()
    for (c in 2:num.cells.per.batch.empty){
        tmp.cell <- cbind(tmp.cell, emptyCell())
    }
    for (c in 1:num.cells.per.batch.bad){
        tmp.cell <- cbind(tmp.cell, badCell())
    }
    bar.bad <- paste(batch.prefix, bar.codes, sep="_")
    colnames(tmp.cell) <- bar.bad
    return(tmp.cell)
}

makeBatch <- function(batch.prefix, odd=TRUE){
    batch.good <- generateGoodBatch(batch.prefix, bar.codes=if(odd) barcodes.1_96 else barcodes.97_192)
    batch.bad <- generateBadBatch(batch.prefix, bar.codes=if(odd) barcodes.97_192 else barcodes.1_96)
    batch <- cbind(batch.good, batch.bad)

    ## Shuffle columns
    batch.shuff <- batch[,sample(ncol(batch))]
    return(batch.shuff)
}

batch.1 <- makeBatch("P1_B1", odd=TRUE)
batch.2 <- makeBatch("P1_B2", odd=FALSE)
batch.3 <- makeBatch("P2_B3", odd=TRUE)
batch.4 <- makeBatch("P2_B4", odd=FALSE)

total.matrix <- cbind(batch.1, batch.2, batch.3, batch.4)
rownames(total.matrix) <- paste("GENE", 1:nrow(total.matrix), sep="")

write.table(as.matrix(total.matrix), file="test.matrix",
            col.names=NA, sep='\t', quote=FALSE)