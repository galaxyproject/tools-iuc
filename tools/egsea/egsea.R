# Code based on (and inspired by) the Galaxy limma-voom/edgeR/DESeq2 wrappers

options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(EGSEA)
    library(limma)
    library(edgeR)
    library(optparse)
})


## Function Declaration

sanitiseEquation <- function(equation) {
    equation <- gsub(" *[+] *", "+", equation)
    equation <- gsub(" *[-] *", "-", equation)
    equation <- gsub(" *[/] *", "/", equation)
    equation <- gsub(" *[*] *", "*", equation)
    equation <- gsub("^\\s+|\\s+$", "", equation)
    return(equation)
}

# Function to sanitise group information
sanitiseGroups <- function(string) {
    string <- gsub(" *[,] *", ",", string)
    string <- gsub("^\\s+|\\s+$", "", string)
    return(string)
}

# Generating design information
pasteListName <- function(string) {
    return(paste0("factors$", string))
}

## Input Processing

option_list <- list(
    make_option(c("-threads", "--threads"), default=2, type="integer", help="Number of threads for egsea"),
    make_option(c("-filesPath", "--filesPath"), type="character", help="JSON list object if multiple files input"),
    make_option(c("-matrixPath", "--matrixPath"), type="character", help="Path to count matrix"),
    make_option(c("-factFile", "--factFile"), type="character", help="Path to factor information file"),
    make_option(c("-factInput", "--factInput"), type="character", help="String containing factors if manually input"),
    make_option(c("-contrastData", "--contrastData"), type="character", help="Contrasts of Interest (Groups to compare)"),
    make_option(c("-genes", "--genes"), type="character", help="Path to genes file"),
    make_option(c("-species", "--species"), type="character"),
    make_option(c("-base_methods", "--base_methods"), type="character", help="Gene set testing methods"),
    make_option(c("-msigdb", "--msigdb"), type="character", help="MSigDB Gene Set Collections"),
    make_option(c("-keggdb", "--keggdb"), type="character", help="KEGG Pathways"),
    make_option(c("-keggupdated", "--keggupdated"), type="logical", help="Use updated KEGG"),
    make_option(c("-gsdb", "--gsdb"), type="character", help = "GeneSetDB Gene Sets"),
    make_option(c("-display_top", "--display_top"), type="integer", help = "Number of top Gene Sets to display"),
    make_option(c("-min_size", "--min_size"), type="integer", help = "Minimum Size of Gene Set"),
    make_option(c("-fdr_cutoff", "--fdr_cutoff"), type="double", help = "FDR cutoff"),
    make_option(c("-combine_method", "--combine_method"), type="character", help="Method to use to combine the p-values"),
    make_option(c("-sort_method", "--sort_method"), type="character", help="Method to sort the results"),
    make_option(c("-rdaOpt", "--rdaOpt"), type="character", help="Output RData file")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)


## Read in Files

if (!is.null(args$filesPath)) {
    # Process the separate count files (adapted from DESeq2 wrapper)
    library("rjson")
    parser <- newJSONParser()
    parser$addData(args$filesPath)
    factorList <- parser$getObject()
    factors <- sapply(factorList, function(x) x[[1]])
    filenamesIn <- unname(unlist(factorList[[1]][[2]]))
    sampleTable <- data.frame(sample=basename(filenamesIn),
                            filename=filenamesIn,
                            row.names=filenamesIn,
                            stringsAsFactors=FALSE)
    for (factor in factorList) {
        factorName <- factor[[1]]
        sampleTable[[factorName]] <- character(nrow(sampleTable))
        lvls <- sapply(factor[[2]], function(x) names(x))
        for (i in seq_along(factor[[2]])) {
            files <- factor[[2]][[i]][[1]]
            sampleTable[files,factorName] <- lvls[i]
        }
        sampleTable[[factorName]] <- factor(sampleTable[[factorName]], levels=lvls)
    }
    rownames(sampleTable) <- sampleTable$sample
    rem <- c("sample","filename")
    factors <- sampleTable[, !(names(sampleTable) %in% rem), drop=FALSE]

    #read in count files and create single table
    countfiles <- lapply(sampleTable$filename, function(x){read.delim(x, row.names=1)})
    counts <- do.call("cbind", countfiles)

} else {
 # Process the single count matrix
    counts <- read.table(args$matrixPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[ , -1]
    countsRows <- nrow(counts)

    # Process factors
    if (is.null(args$factInput)) {
            factorData <- read.table(opt$factFile, header=TRUE, sep="\t", strip.white=TRUE)
            # order samples as in counts matrix
            factorData <- factorData[match(colnames(counts), factorData[, 1]), ]
            factors <- factorData[, -1, drop=FALSE]
    }  else {
            factors <- unlist(strsplit(args$factInput, "|", fixed=TRUE))
            factorData <- list()
            for (fact in factors) {
                newFact <- unlist(strsplit(fact, split="::"))
                factorData <- rbind(factorData, newFact)
            } # Factors have the form: FACT_NAME::LEVEL,LEVEL,LEVEL,LEVEL,... The first factor is the Primary Factor.

            # Set the row names to be the name of the factor and delete first row
            row.names(factorData) <- factorData[, 1]
            factorData <- factorData[, -1]
            factorData <- sapply(factorData, sanitiseGroups)
            factorData <- sapply(factorData, strsplit, split=",")
            factorData <- sapply(factorData, make.names)
            # Transform factor data into data frame of R factor objects
            factors <- data.frame(factorData)
    }
}

# Create a DGEList object
counts <- DGEList(counts)

# Set group to be the Primary Factor input
group <- factors[, 1, drop=FALSE]

# Split up contrasts separated by comma into a vector then sanitise
contrastData <- unlist(strsplit(args$contrastData, split=","))
contrastData <- sanitiseEquation(contrastData)
contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)

# Creating design
row.names(factors) <- colnames(counts)
factorList <- sapply(names(factors), pasteListName)

formula <- "~0"
for (i in 1:length(factorList)) {
    formula <- paste(formula, factorList[i], sep="+")
}
formula <- formula(formula)

design <- model.matrix(formula)

for (i in 1:length(factorList)) {
    colnames(design) <- gsub(factorList[i], "", colnames(design), fixed=TRUE)
}

## Generate Contrasts information
contrasts <- makeContrasts(contrasts=contrastData, levels=design)


## Add Gene Symbol information

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

gs.annots <- buildIdx(entrezIDs=rownames(counts), species=args$species, msigdb.gsets=msigdb, gsdb.gsets=gsdb, kegg.exclude=kegg_exclude, kegg.updated=args$keggupdated)


## Run egsea.cnt

gsa <- egsea.cnt(counts=counts, group=group, design=design, contrasts=contrasts, gs.annots=gs.annots, symbolsMap=genes, baseGSEAs=base_methods, minSize=args$min_size, display.top=args$display_top, combineMethod=args$combine_method, sort.by=args$sort_method, report.dir='./report_dir', fdr.cutoff=args$fdr_cutoff, num.threads=args$threads, report=TRUE)


## Output RData file

if (!is.null(args$rdaOpt)) {
  save.image(file = "EGSEA_analysis.RData")
}