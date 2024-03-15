#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list(
    make_option(c("--BIOMfilename"), action = "store", dest = "biom", help = "Input BIOM file"),
    make_option(c("--treefilename"), action = "store", dest = "tree", default = NULL, help = "Input Tree newick/nexus file"),
    make_option(c("--parseFunction"), action = "store", dest = "parsefoo", default = "parse_taxonomy_default", help = "Parse function parse_taxonomy_default/read_tree_greengenes"),
    make_option(c("--refseqfilename"), action = "store", dest = "sequences", default = NULL, help = "Input Sequence fasta file"),
    make_option(c("--output"), action = "store", dest = "output", help = "RDS output")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

parsefoo <- get(opt$parsefoo)
phyloseq_obj <- import_biom(
    BIOMfilename = opt$biom,
    treefilename = opt$tree,
    refseqfilename = opt$sequences,
    parseFunction=parsefoo
)

print(phyloseq_obj)

# save R object to file
saveRDS(phyloseq_obj, file = opt$output, compress = TRUE)
