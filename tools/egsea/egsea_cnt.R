# Code based on (and inspired by) the Galaxy limma-voom/edgeR/DESeq2 wrappers

options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr()); q("no", 1, F)
  })

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(EGSEA)
    library(limma)
    library(edgeR)
    library(optparse)
})


# Function Declaration
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
paste_list_name <- function(string) {
    return(paste0("factors$", string))
}

## Input Processing

option_list <- list(
    make_option(c("-t", "--threads"),
                type = "integer",
                default = 2,
                help = "Number of threads for egsea",
                metavar = "integer"),
    make_option(c("-f", "--filesPath"),
                type = "character",
                help = "JSON list object if multiple files input",
                metavar = "character"),
    make_option(c("-m", "--matrixPath"),
                type = "character",
                help = "Path to count matrix",
                metavar = "character"),
    make_option(c("-i", "--factFile"),
                type = "character",
                help = "Path to factor information file",
                metavar = "character"),
    make_option(c("-x", "--factInput"),
                type = "character",
                help = "String containing factors if manually input",
                metavar = "character"),
    make_option(c("-c", "--contrast_data"),
                type = "character",
                help = "Contrasts of Interest (Groups to compare)",
                metavar = "character"),
    make_option(c("-g", "--genes"),
                type = "character",
                help = "Path to genes file",
                metavar = "character"),
    make_option(c("-h", "--species"),
                type = "character",
                metavar = "character"),
    make_option(c("-b", "--base_methods"),
                type = "character",
                help = "Gene set testing methods",
                metavar = "character"),
    make_option(c("-z", "--msigdb"),
                type = "character",
                help = "MSigDB Gene Set Collections",
                metavar = "character"),
    make_option(c("-k", "--keggdb"),
                type = "character",
                help = "KEGG Pathways",
                metavar = "character"),
    make_option(c("-u", "--keggupdated"),
                type = "logical",
                help = "Use updated KEGG",
                metavar = "logical"),
    make_option(c("-y", "--gsdb"),
                type = "character",
                help = "GeneSetDB Gene Sets",
                metavar = "character"),
    make_option(c("-t", "--display_top"),
                type = "integer",
                help = "Number of top Gene Sets to display",
                metavar = "character"),
    make_option(c("-l", "--min_size"),
                type = "integer",
                help = "Minimum Size of Gene Set",
                metavar = "character"),
    make_option(c("-z", "--fdr_cutoff"),
                type = "double",
                help = "FDR cutoff",
                metavar = "character"),
    make_option(c("-p", "--combine_method"),
                type = "character",
                help = "Method to use to combine the p-values",
                metavar = "character"),
    make_option(c("-s", "--sort_method"),
                type = "character",
                help = "Method to sort the results",
                metavar = "character"),
    make_option(c("-r", "--rdaOpt"),
                type = "character",
                help = "Output RData file",
                metavar = "character")
    )

opt_parser <- OptionParser(usage = "%prog [options] file",
                           option_list = option_list)
opt <- parse_opt(opt_parser)


## Read in Files

if (!is.null(opt$filesPath)) {
    # Process the separate count files (adapted from DESeq2 wrapper)
    library("rjson")
    parser <- newJSONParser()
    parser$addData(opt$filesPath)
    factor_list <- parser$getObject()
    factors <- sapply(factor_list, function(x) x[[1]])
    filenames_in <- unname(unlist(factor_list[[1]][[2]]))
    sample_table <- data.frame(sample = basename(filenames_in),
                            filename = filenames_in,
                            row.names = filenames_in,
                            stringsAsFactors = FALSE)
    for (factor in factor_list) {
        factor_name <- factor[[1]]
        sample_table[[factor_name]] <- character(nrow(sample_table))
        lvls <- sapply(factor[[2]], function(x) names(x))
        for (i in seq_along(factor[[2]])) {
            files <- factor[[2]][[i]][[1]]
            sample_table[files, factor_name] <- lvls[i]
        }
        sample_table[[factor_name]] <- factor(sample_table[[factor_name]],
                                              levels = lvls)
    }
    rownames(sample_table) <- sample_table$sample
    rem <- c("sample", "filename")
    factors <- sample_table[, !(names(sample_table) %in% rem), drop = FALSE]

    #read in count files and create single table
    countfiles <- lapply(sample_table$filename, function(x) {
      read.delim(x, row.names = 1)
      })
    counts <- do.call("cbind", countfiles)

} else {
 # Process the single count matrix
    counts <- read.table(opt$matrixPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[, -1]
    counts_rows <- nrow(counts)

    # Process factors
    if (is.null(opt$factInput)) {
            factor_data <- read.table(opt$factFile, header = TRUE, sep = "\t",
                                     strip.white = TRUE)
            # check samples names match
            if (!any(factor_data[, 1] %in% colnames(counts))) {
                stop("Sample IDs in factors file and count matrix don't match")
            }
            # order samples as in counts matrix
            factor_data <- factor_data[match(colnames(counts),
                                             factor_data[, 1]), ]
            factors <- factor_data[, -1, drop = FALSE]
    }  else {
            factors <- unlist(strsplit(opt$factInput, "|", fixed = TRUE))
            factor_data <- list()
            for (fact in factors) {
                new_fact <- unlist(strsplit(fact, split = "::"))
                factor_data <- rbind(factor_data, new_fact)
            }
            # Factors have the form: FACT_NAME::LEVEL,LEVEL,LEVEL,LEVEL,...
            # The first factor is the Primary Factor.

            # Set the row names to be the name of the factor and
            # delete first row
            row.names(factor_data) <- factor_data[, 1]
            factor_data <- factor_data[, -1]
            factor_data <- sapply(factor_data, sanitise_groups)
            factor_data <- sapply(factor_data, strsplit, split = ",")
            factor_data <- sapply(factor_data, make.names)
            # Transform factor data into data frame of R factor objects
            factors <- data.frame(factor_data)
    }
}

# Create a DGEList object
counts <- DGEList(counts)

# Set group to be the Primary Factor input
group <- factors[, 1, drop = FALSE]

# Split up contrasts separated by comma into a vector then sanitise
contrast_data <- unlist(strsplit(opt$contrast_data, split = ","))
contrast_data <- sanitise_equation(contrast_data)
contrast_data <- gsub(" ", ".", contrast_data, fixed = TRUE)

# Creating design
row.names(factors) <- colnames(counts)
factor_list <- sapply(names(factors), paste_list_name)

formula <- "~0"
for (i in seq_len(factor_list)) {
    formula <- paste(formula, factor_list[i], sep = "+")
}
formula <- formula(formula)

design <- model.matrix(formula)

for (i in seq_len(factor_list)) {
    colnames(design) <- gsub(factor_list[i], "", colnames(design), fixed = TRUE)
}

## Generate Contrasts information
contrasts <- makeContrasts(contrasts = contrast_data, levels = design)


## Add Gene Symbol information

genes <- read.table(opt$genes, sep = "\t", header = TRUE)


## Set Gene Set Testing Methods

base_methods <- unlist(strsplit(opt$base_methods, ","))


## Set Gene Sets

if (opt$msigdb != "None") {
    msigdb <- unlist(strsplit(opt$msigdb, ","))
} else {
    msigdb <- "none"
}

if (opt$keggdb != "None") {
    keggdb <- unlist(strsplit(opt$keggdb, ","))
    kegg_all <- c("Metabolism" = "keggmet",
                  "Signaling" = "keggsig",
                  "Disease" = "keggdis")
    kegg_exclude <- names(kegg_all[!(kegg_all %in% keggdb)])
} else {
    kegg_exclude <- "all"
}

if (opt$gsdb != "None") {
    gsdb <- unlist(strsplit(opt$gsdb, ","))
} else {
    gsdb <- "none"
}

## Index gene sets

gs_annots <- buildIdx(entrezIDs = rownames(counts),
                      species = opt$species,
                      msigdb.gsets = msigdb,
                      gsdb.gsets = gsdb,
                      kegg.exclude = kegg_exclude,
                      kegg.updated = opt$keggupdated)


## Run egsea.cnt

gsa <- egsea.cnt(counts = counts,
                 group = group,
                 design = design,
                 contrasts = contrasts,
                 gs_annots = gs_annots,
                 symbolsMap = genes,
                 baseGSEAs = base_methods,
                 minSize = opt$min_size,
                 display.top = opt$display_top,
                 combineMethod = opt$combine_method,
                 sort.by = opt$sort_method,
                 report.dir = "./report_dir",
                 fdr.cutoff = opt$fdr_cutoff,
                 num.threads = opt$threads,
                 report = TRUE)


## Output RData file

if (!is.null(opt$rdaOpt)) {
  save.image(file = "EGSEA_analysis.RData")
}
