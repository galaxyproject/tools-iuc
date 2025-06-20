get_deseq_dataset <- function(sample_table, header, design_formula, tximport, txtype, tx2gene) {
    dir <- ""

    has_header <- !is.null(header)
    use_txi <- !is.null(tximport)
    if (use_txi) {
        if (is.null(tx2gene)) stop("A transcript-to-gene map or a GTF/GFF3 file is required for tximport")
        if (tolower(file_ext(tx2gene)) == "gff") {
            gff_file <- tx2gene
        } else {
            gff_file <- NULL
            tx2gene <- read.table(tx2gene, header = has_header)
        }
    }

    if (!use_txi && has_header) {
        countfiles <- lapply(as.character(sample_table$filename), read.delim, row.names = 1)
        tbl <- do.call("cbind", countfiles)
        colnames(tbl) <- rownames(sample_table) # take sample ids from header

        # check for htseq report lines (from DESeqDataSetFromHTSeqCount function)
        old_special_names <- c(
            "no_feature",
            "ambiguous",
            "too_low_aQual",
            "not_aligned",
            "alignment_not_unique"
        )
        special_rows <- (substr(rownames(tbl), 1, 1) == "_") | rownames(tbl) %in% old_special_names
        tbl <- tbl[!special_rows, , drop = FALSE]

        dds <- DESeqDataSetFromMatrix(
            countData = tbl,
            colData = subset(sample_table, select = -filename),
            design = design_formula
        )
    } else if (!use_txi && !has_header) {
        # construct the object from HTSeq files
        dds <- DESeqDataSetFromHTSeqCount(
            sampleTable = sample_table,
            directory = dir,
            design = design_formula
        )
        colnames(dds) <- row.names(sample_table)
    } else {
        # construct the object using tximport
        library("tximport")
        txi_files <- as.character(sample_table$filename)
        labs <- row.names(sample_table)
        names(txi_files) <- labs
        if (!is.null(gff_file)) {
            # first need to make the tx2gene table
            # this takes ~2-3 minutes using Bioconductor functions
            suppressPackageStartupMessages({
                library("GenomicFeatures")
            })
            txdb <- makeTxDbFromGFF(gff_file)
            k <- keys(txdb, keytype = "TXNAME")
            tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
            # Remove 'transcript:' from transcript IDs (when gff_file is a GFF3 from Ensembl and the transcript does not have a Name)
            tx2gene$TXNAME <- sub("^transcript:", "", tx2gene$TXNAME) # nolint
        }
        try(txi <- tximport(txi_files, type = txtype, tx2gene = tx2gene))
        if (!exists("txi")) {
            # Remove version from transcript IDs in tx2gene...
            tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME) # nolint
            # ...and in txi_files
            txi <- tximport(txi_files, type = txtype, tx2gene = tx2gene, ignoreTxVersion = TRUE)
        }
        dds <- DESeqDataSetFromTximport(
            txi,
            subset(sample_table, select = -c(filename)),
            design_formula
        )
    }
    return(dds)
}
