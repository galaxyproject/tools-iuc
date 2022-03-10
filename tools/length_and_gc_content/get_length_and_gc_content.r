# originally by Devon Ryan, https://www.biostars.org/p/84467/

options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("GenomicRanges")
    library("rtracklayer")
    library("Rsamtools")
    library("optparse")
    library("data.table")
})

option_list <- list(
    make_option(c("-g","--gtf"), type="character", help="Input GTF file with gene / exon information."),
    make_option(c("-f","--fasta"), type="character", default=FALSE, help="FASTA file that corresponds to the supplied GTF."),
    make_option(c("-l","--length"), type="character", default=FALSE, help="Output file with Gene ID and length."),
    make_option(c("-c","--gc_content"), type="character", default=FALSE, help="Output file with Gene ID and GC content.")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

GTFfile = args$gtf
FASTAfile = args$fasta
length = args$length
gc_content = args$gc_content

# Check args:
if (is.logical(FASTAfile) & ! is.logical(gc_content)){
    stop("gc_content output requires FASTA input")
}
if (is.logical(length) & is.logical(gc_content)){
    stop("neither gc_content nor length was set nothing to do.")
}

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome=NA, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

if (! is.logical(FASTAfile)){
    #Open the fasta file
    FASTA <- FaFile(FASTAfile)
    open(FASTA)

    #Add the GC numbers
    elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
}
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
if (! is.logical(gc_content)){
    calc_GC_length <- function(x) {
        nGCs = sum(elementMetadata(x)$nGCs)
        width = sum(elementMetadata(x)$widths)
        c(width, nGCs/width)
    }
    output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
    output <- data.frame(setDT(data.frame(output), keep.rownames = TRUE)[])
} else {
    all_widths <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), function(x){sum(elementMetadata(x)$widths)})
    output <- data.frame(gene_id = names(all_widths), length = all_widths, gc_content = NA)
}


write.table(output[,c(1,2)], file=length, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(output[,c(1,3)], file=gc_content, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


sessionInfo()
