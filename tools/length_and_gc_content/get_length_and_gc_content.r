# originally by Devon Ryan, https://www.biostars.org/p/84467/

options(show.error.messages = F,
        error = function() {
          cat(geterrmessage(), file = stderr())
          q("no", 1, F)
        })

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
  make_option(c("-g", "--gtf"), type = "character",
              help = "Input gtf file with gene / exon information."),
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "fasta file that corresponds to the supplied gtf."),
  make_option(c("-l", "--length"), type = "character", default = NULL,
              help = "Output file with Gene ID and length."),
  make_option(c("-c", "--gc_content"), type = "character", default = NULL,
              help = "Output file with Gene ID and GC content.")
)

parser <- OptionParser(usage = "%prog [options] file",
                       option_list = option_list)
args <- parse_args(parser)

gtf_file <- args$gtf
fasta_file <- args$fasta
length <- args$length
gc_content <- args$gc_content

# Check args:
if (is.null(fasta_file) & !is.null(gc_content)) {
  stop("gc_content output requires fasta input")
}
if (is.null(length) & is.null(gc_content)) {
  stop("neither gc_content nor length was set nothing to do.")
}

#Load the annotation and reduce it
gtf <- import.gff(gtf_file, format = "gtf", genome = NA, feature.type = "exon")
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_id))
reduced_gtf <- unlist(grl, use.names = T)
elementMetadata(reduced_gtf)$gene_id <- rep(names(grl), elementNROWS(grl))

if (! is.null(gc_content)) {
  #Open the fasta file
  fasta <- FaFile(fasta_file)
  open(fasta)

  #Add the GC numbers
  elementMetadata(reduced_gtf)$n_gcs <-
    letterFrequency(getSeq(fasta, reduced_gtf), "GC")[, 1]
}
elementMetadata(reduced_gtf)$widths <- width(reduced_gtf)

#Create a list of the ensembl_id/GC/length
if (! is.null(gc_content)) {
  calc_gc_length <- function(x) {
    n_gcs <- sum(elementMetadata(x)$n_gcs)
    width <- sum(elementMetadata(x)$widths)
    c(width, n_gcs / width)
  }
  output <- t(sapply(split(reduced_gtf, elementMetadata(reduced_gtf)$gene_id),
                     calc_gc_length))
  output <- data.frame(setDT(data.frame(output), keep.rownames = TRUE)[])
  write.table(output[, c(1, 3)], file = gc_content,
              col.names = FALSE, row.names = FALSE,
              quote = FALSE, sep = "\t")
} else {
  all_widths <- sapply(split(reduced_gtf, elementMetadata(reduced_gtf)$gene_id),
                       function(x) {
                         sum(elementMetadata(x)$widths)
                        })
  output <- data.frame(gene_id = names(all_widths),
                       length = all_widths)
}

if (! is.null(length)) {
  write.table(output[, c(1, 2)], file = length,
              col.names = FALSE, row.names = FALSE,
              quote = FALSE, sep = "\t")
}


sessionInfo()
