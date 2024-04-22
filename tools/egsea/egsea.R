# Code based on (and inspired by) the Galaxy limma-voom/edgeR/DESeq2 wrappers

options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(EGSEA)
    library(limma)
    library(edgeR)
    library(optparse)
})


## Function Declaration

sanitise_equation <- function(equation) {
    equation <- gsub(" *[+] *", "+", equation)
    equation <- gsub(" *[-] *", "-", equation)
    equation <- gsub(" *[/] *", "/", equation)
    equation <- gsub(" *[*] *", "*", equation)
    equation <- gsub("^\\s+|\\s+$", "", equation)
    return(equation)
}

# Function to sanitise group information
sanitise_groups <- function(string) {
    string <- gsub(" *[,] *", ",", string)
    string <- gsub("^\\s+|\\s+$", "", string)
    return(string)
}

# Generating design information
paste_listname <- function(string) {
    return(paste0("factors$", string))
}

## Input Processing

option_list <- list(
    make_option("--threads", default = 2, type = "integer", help = "Number of threads for egsea"),
    make_option("--filesPath", type = "character", help = "JSON list object if multiple files input"),
    make_option("--matrixPath", type = "character", help = "Path to count matrix"),
    make_option("--factFile", type = "character", help = "Path to factor information file"),
    make_option("--factInput", type = "character", help = "String containing factors if manually input"),
    make_option("--contrastData", type = "character", help = "Contrasts of Interest (Groups to compare)"),
    make_option("--genes", type = "character", help = "Path to genes file"),
    make_option("--species", type = "character"),
    make_option("--base_methods", type = "character", help = "Gene set testing methods"),
    make_option("--msigdb", type = "character", help = "MSigDB Gene Set Collections"),
    make_option("--keggdb", type = "character", help = "KEGG Pathways"),
    make_option("--keggupdated", type = "logical", help = "Use updated KEGG"),
    make_option("--gsdb", type = "character", help = "GeneSetDB Gene Sets"),
    make_option("--display_top", type = "integer", help = "Number of top Gene Sets to display"),
    make_option("--min_size", type = "integer", help = "Minimum Size of Gene Set"),
    make_option("--fdr_cutoff", type = "double", help = "FDR cutoff"),
    make_option("--combine_method", type = "character", help = "Method to use to combine the p-values"),
    make_option("--sort_method", type = "character", help = "Method to sort the results"),
    make_option("--rdaOpt", type = "character", help = "Output RData file")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


## Read in Files

if (!is.null(args$filesPath)) {
    # Process the separate count files (adapted from DESeq2 wrapper)
    library("rjson")
    parser <- newJSONParser()
    parser$addData(args$filesPath)
    factor_list <- parser$getObject()
    factors <- sapply(factor_list, function(x) x[[1]])
    filenames_in <- unname(unlist(factor_list[[1]][[2]]))
    sampletable <- data.frame(
        sample = basename(filenames_in),
        filename = filenames_in,
        row.names = filenames_in,
        stringsAsFactors = FALSE
    )
    for (factor in factor_list) {
        factorname <- factor[[1]]
        sampletable[[factorname]] <- character(nrow(sampletable))
        lvls <- sapply(factor[[2]], function(x) names(x))
        for (i in seq_along(factor[[2]])) {
            files <- factor[[2]][[i]][[1]]
            sampletable[files, factorname] <- lvls[i]
        }
        sampletable[[factorname]] <- factor(sampletable[[factorname]], levels = lvls)
    }
    rownames(sampletable) <- sampletable$sample
    rem <- c("sample", "filename")
    factors <- sampletable[, !(names(sampletable) %in% rem), drop = FALSE]

    # read in count files and create single table
    countfiles <- lapply(sampletable$filename, function(x) {
        read.delim(x, row.names = 1)
    })
    counts <- do.call("cbind", countfiles)
} else {
    # Process the single count matrix
    counts <- read.table(args$matrixPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[, -1]
    countsrows <- nrow(counts)

    # Process factors
    if (is.null(args$factInput)) {
        factordata <- read.table(args$factFile, header = TRUE, sep = "\t", strip.white = TRUE, stringsAsFactors = TRUE)
        # check samples names match
        if (!any(factordata[, 1] %in% colnames(counts))) {
            stop("Sample IDs in factors file and count matrix don't match")
        }
        # order samples as in counts matrix
        factordata <- factordata[match(colnames(counts), factordata[, 1]), ]
        factors <- factordata[, -1, drop = FALSE]
    } else {
        factors <- unlist(strsplit(args$factInput, "|", fixed = TRUE))
        factordata <- list()
        for (fact in factors) {
            newfact <- unlist(strsplit(fact, split = "::"))
            factordata <- rbind(factordata, newfact)
        } # Factors have the form: FACT_NAME::LEVEL,LEVEL,LEVEL,LEVEL,... The first factor is the Primary Factor.

        # Set the row names to be the name of the factor and delete first row
        row.names(factordata) <- factordata[, 1]
        factordata <- factordata[, -1]
        factordata <- sapply(factordata, sanitise_groups)
        factordata <- sapply(factordata, strsplit, split = ",")
        factordata <- sapply(factordata, make.names)
        # Transform factor data into data frame of R factor objects
        factors <- data.frame(factordata, stringsAsFactors = TRUE)
    }
}

# Create a DGEList object
counts <- DGEList(counts)

# Set group to be the Primary Factor input
group <- factors[, 1, drop = FALSE]

# Split up contrasts separated by comma into a vector then sanitise
contrast_data <- unlist(strsplit(args$contrastData, split = ","))
contrast_data <- sanitise_equation(contrast_data)
contrast_data <- gsub(" ", ".", contrast_data, fixed = TRUE)

# Creating design
row.names(factors) <- colnames(counts)
factor_list <- sapply(names(factors), paste_listname)

formula <- "~0"
for (i in seq_along(factor_list)) {
    formula <- paste(formula, factor_list[i], sep = "+")
}
formula <- formula(formula)

design <- model.matrix(formula)

for (i in seq_along(factor_list)) {
    colnames(design) <- gsub(factor_list[i], "", colnames(design), fixed = TRUE)
}

## Generate Contrasts information
contrasts <- makeContrasts(contrasts = contrast_data, levels = design)


## Add Gene Symbol information

genes <- read.table(args$genes, sep = "\t", header = TRUE)


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
    kegg_all <- c("Metabolism" = "keggmet", "Signaling" = "keggsig", "Disease" = "keggdis")
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

gs_annots <- buildIdx(entrezIDs = rownames(counts), species = args$species, msigdb.gsets = msigdb, gsdb.gsets = gsdb, kegg.exclude = kegg_exclude, kegg.updated = args$keggupdated)


## Run egsea.cnt

gsa <- egsea.cnt(counts = counts, group = group, design = design, contrasts = contrasts, gs.annots = gs_annots, symbolsMap = genes, baseGSEAs = base_methods, minSize = args$min_size, display.top = args$display_top, combineMethod = args$combine_method, sort.by = args$sort_method, report.dir = "./report_dir", fdr.cutoff = args$fdr_cutoff, num.threads = args$threads, report = TRUE)


## Output RData file

if (!is.null(args$rdaOpt)) {
    save.image(file = "EGSEA_analysis.RData")
}
