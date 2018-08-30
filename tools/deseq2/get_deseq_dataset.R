get_deseq_dataset <- function(sampleTable, header, designFormula, tximport, txtype, tx2gene) {

  dir <- ""

  if (!is.null(header)) {
    hasHeader <- TRUE
  } else {
    hasHeader <- FALSE
  }

  if (!is.null(tximport)) {
    if (is.null(tx2gene)) stop("A transcript-to-gene map or a GTF file is required for tximport")
    if (tolower(file_ext(opt$tx2gene)) == "gtf") {
      gtfFile <-tx2gene
    } else {
      gtfFile <- NULL
      tx2gene <- read.table(tx2gene, header=FALSE)
    }
    useTXI <- TRUE
  } else {
    useTXI <- FALSE
  }

  if (!useTXI & hasHeader) {
      countfiles <- lapply(as.character(sampleTable$filename), function(x){read.delim(x, row.names=1)})
      tbl <- do.call("cbind", countfiles)
      colnames(tbl) <- rownames(sampleTable) # take sample ids from header

      # check for htseq report lines (from DESeqDataSetFromHTSeqCount function)
      oldSpecialNames <- c("no_feature", "ambiguous", "too_low_aQual",
          "not_aligned", "alignment_not_unique")
      specialRows <- (substr(rownames(tbl), 1, 1) == "_") | rownames(tbl) %in% oldSpecialNames
      tbl <- tbl[!specialRows, , drop = FALSE]

      dds <- DESeqDataSetFromMatrix(countData = tbl,
                                    colData = subset(sampleTable, select=-(filename)),
                                    design = designFormula)
  } else if (!useTXI & !hasHeader) {

    # construct the object from HTSeq files
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dir,
                                      design =  designFormula)
    colnames(dds) <- row.names(sampleTable)

  } else {
      # construct the object using tximport
      # first need to make the tx2gene table
      # this takes ~2-3 minutes using Bioconductor functions
      if (!is.null(gtfFile)) {
        suppressPackageStartupMessages({
          library("GenomicFeatures")
        })
        txdb <- makeTxDbFromGFF(gtfFile, format="gtf")
        k <- keys(txdb, keytype = "GENEID")
        df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
        tx2gene <- df[, 2:1]  # tx ID, then gene ID
      }
      library("tximport")
      txiFiles <- as.character(sampleTable$filename)
      labs <- row.names(sampleTable)
      names(txiFiles) <- labs
      txi <- tximport(txiFiles, type=txtype, tx2gene=tx2gene)
      dds <- DESeqDataSetFromTximport(txi,
                                      subset(sampleTable, select=-c(filename)),
                                      designFormula)
  }
  return(dds)
}
