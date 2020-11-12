options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr()); q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  
  library(EGSEA)
  library(limma)
  library(edgeR)
  
})

option_list <- list(
  make_option(c("-g", "--geneIDs"),
              type = "character",
              help = "Character, a vector of Gene IDs to be tested for ORA.
              They must be Entrez IDs if EGSEAdata collections are used.",
              metavar = "character"),
  make_option(c("-u", "--universe"),
              type = "Character",
              help = "character, a vector of Enterz IDs to be used as a background list.
              If universe=NULL, the background list is created
              from the AnnotationDbi package.",
              metavar = "Character"),
  make_option(c("-l", "--logFC"),
              type = "double",
              help = "double,  it can be a matrix or vector of the same length of entrezIDs. If logFC=NULL, 1 is used as a default value. Then, the regulation direction in heatmaps and pathway maps is not indicative of the gene regulation direction.",
              metavar = "double"),
  make_option(c("-t", "--titel"),
              type = "character",
              help = "character, a short description of the experimental contrast.",
              metavar = "character"),
  make_option(c("-a", "--gs.annots"),
              type = "list",
              help = "list, list of objects of class GSCollectionIndex. It is generated using one of these functions: buildIdx, buildMSigDBIdx, buildKEGGIdx, buildGeneSetDBIdx, and buildCustomIdx.",
              metavar = "list"),
  make_option(c("-s", "--symbolsMap"),
              type = "dataframe",
              help = "dataframe, an K x 2 matrix stores the gene symbol of each Entrez Gene ID. The first column must be the Entrez Gene IDs and the second column must be the Gene Symbols. It is used for the heatmap visualization.",
              metavar = "dataframe"),
  make_option(c("-m", "--minSize"),
              type = "integer",
              help = "integer, the minimum size of a gene set to be included in the analysis. Default minSize= 2.",
              metavar = "integer"),
  make_option(c("-d", "--display.top"),
              type = "integer",
              help = "integer, the number of top gene sets to be displayed in the EGSEA report.",
              metavar = "integer"),
  make_option(c("-x", "--sort.by"),
              type = "character",
              help = "character, determines how to order the analysis results in the stats table.",
              metavar = "character"),
  make_option(c("-r", "--report.dir"),
              type = "character",
              help = "character, directory into which the analysis results are written out.",
              metavar = "character"),
  make_option(c("-k", "--kegg.dir"),
              type = "character",
              help = "character, the directory of KEGG pathway data file (.xml) and image file (.png).",
              metavar = "character"),
  make_option(c("-p", "--sum.plot.axis"),
              type = "character",
              help = "character, the x-axis of the summary plot. All the values accepted by the sort.by parameter can be used.",
              metavar = "character"),
  make_option(c("-c", "--sum.plot.cutoff"),
              type = "integer",
              help = "numeric, cut-off threshold to filter the gene sets of the summary plots based on the values of the sum.plot.axis.",
              metavar = "integer"),
  make_option(c("-n", "--num.threads"),
              type = "integer",
              help = "numeric, number of CPU cores to be used.",
              metavar = "integer"),
  make_option(c("-y", "--report"),
              type = "logical",
              help = "logical, whether to generate the EGSEA interactive report. It takes longer time to run. Default is True.",
              metavar = "logical"),
  make_option(c("-i", "--interactive"),
              type = "logical",
              help = "logical, whether to generate interactive tables and plots. Note this might dramatically increase the size of the EGSEA report.",
              metavar = "logical"),
  make_option(c("-v", "--verbose"),
              type = "logical",
              help = "logical, whether to print out progress messages and warnings.",
              metavar = "logical")
)

opt_parser <- OptionParser(usage = "%prog [options] file",
                           option_list = option_list)
opt <- parse_args(opt_parser)


# Create a vector of Gene IDs to be tested for ORA
geneIDs <- grep("^#", readLines(opt$geneIDs), invert = T, value = T)

# Create a vector of Enterz IDs to be used as a background list
universe <- c()
for (i in opt$universe) {
  pathways_out <- fromJSON(paste(readLines(i), collapse = ""))
  pathways_list <- pathways_out[[2]]
  names(pathways_list) <- unlist(pathways_out[[1]])
  universe <- unique(c(unlist(pathways_list), universe))
}

# 

logFC <- grep("^#", readLines(opt$logFC), invert = T, value = T)




egsea.ora(geneIDs = geneIDs, universe = universe, logFC = NULL, title = NULL, gs.annots,
          symbolsMap = NULL, minSize = 2, display.top = 20, sort.by = "p.adj",
          report.dir = NULL, kegg.dir = NULL, sum.plot.axis = "p.adj",
          sum.plot.cutoff = NULL, num.threads = 4, report = TRUE,
          interactive = FALSE, verbose = FALSE)
