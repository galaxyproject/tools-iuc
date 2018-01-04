options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(EGSEA)
    library(optparse)
})


option_list <- list(
    make_option(c("-counts", "--counts"), type="character", help="Path to counts file"),
    make_option(c("-group", "--group"), type="character", help="Path to groups file"),
    make_option(c("-genes", "--genes"), type="character", help="Path to genes file"),
    make_option(c("-species", "--species"), type="character"),
    make_option(c("-base_methods", "--base_methods"), type="character", help="Gene set testing methods"),
    make_option(c("-msigdb", "--msigdb"), type="character", help="MSigDB Gene Set Collections"),
    make_option(c("-keggdb", "--keggdb"), type="character", help="KEGG Pathways"),
    make_option(c("-gsdb", "--gsdb"), type="character", help = "GeneSetDB Gene Sets"),
    make_option(c("-display_top", "--display_top"), type="integer", help = "Number of top Gene Sets to display"),
    make_option(c("-min_size", "--min_size"), type="integer", help = "Minimum Size of Gene Set"),
    make_option(c("-fdr_cutoff", "--fdr_cutoff"), type="double", help = "FDR cutoff"),
    make_option(c("-combine_method", "--combine_method"), type="character", help="Method to use to combine the p-values"),
    make_option(c("-sort_method", "--sort_method"), type="character", help="Method to sort the results"),
    make_option(c("-rdata", "--rdaOpt"), type="character", help="Output RData file")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)


## Read in Files

counts <- read.table(args$counts, sep='\t', row.names=1, header=TRUE)
counts <- as.matrix(counts)
group <- read.table(args$group, sep='\t')
group <- factor(group[,1])
genes <- read.table(args$genes, sep='\t', header=TRUE)


## Set Gene Set Testing Methods

base_methods <- unlist(strsplit(args$base_methods, ","))


## Set Gene Sets

if (args$msigdb != "None") {
    msigdb <- unlist(strsplit(args$msigdb, ","))
} else {
    msigdb <- "none"
}

if (args$keggdb != "None") {
    keggdb <- unlist(strsplit(args$keggdb, ","))
    kegg_all <- c("Metabolism"="keggmet", "Signaling"="keggsig", "Disease"="keggdis")
    kegg_exclude <- names(kegg_all[!(kegg_all %in% keggdb)])
} else {
    kegg_exclude <- "all"
}

if (args$gsdb != "None") {
    gsdb <- unlist(strsplit(args$gsdb, ","))
} else {
    gsdb <- "none"
}


## Index gene sets

gs.annots <- buildIdx(entrezIDs=rownames(counts), species=args$species, msigdb.gsets=msigdb, gsdb.gsets=gsdb, kegg.exclude=kegg_exclude)


## Run egsea.cnt

gsa <- egsea.cnt(counts=counts, group=group, gs.annots=gs.annots, symbolsMap=genes, baseGSEAs=base_methods, minSize=args$min_size, display.top=args$display_top, combineMethod=args$combine_method, sort.by=args$sort_method, report.dir='./report_dir', fdr.cutoff=args$fdr_cutoff, num.threads=2, report=TRUE)


# Output RData file

if (!is.null(args$rdata)) {
  save.image(file = "EGSEA_analysis.RData")
}