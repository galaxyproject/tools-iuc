library(multiGSEA,
    quietly = TRUE,
    warn.conflicts = FALSE
)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)

################################################################################
### Input Processing
################################################################################


# Collect arguments from command line
parser <- ArgumentParser(description = "multiGSEA R script")

parser$add_argument("--transcriptomics",
    required = FALSE,
    help = "Transcriptomics data"
)
parser$add_argument(
    "--transcriptome_ids",
    required = FALSE,
    help = "Transcriptomics ids",
    default = "SYMBOL"
)
parser$add_argument("--proteomics",
    required = FALSE,
    help = "Proteomics data"
)
parser$add_argument(
    "--proteome_ids",
    required = FALSE,
    help = "Proteomics ids",
    default = "SYMBOL"
)
parser$add_argument("--metabolomics",
    required = FALSE,
    help = "Metabolomics data"
)
parser$add_argument(
    "--metabolome_ids",
    required = FALSE,
    help = "Metabolomics ids",
    default = "HMDB"
)
parser$add_argument("--organism",
    required = TRUE,
    help = "Organism"
)
parser$add_argument("--combine_pvalues",
    required = TRUE,
    help = "Combine p-values method"
)
parser$add_argument("--padj_method",
    required = TRUE,
    help = "P-adjustment method"
)
parser$add_argument("--databases",
    required = TRUE,
    help = "Pathway databases"
)

args <- parser$parse_args()

## ----Load library-------------------------------------------------------------

organism_mapping <- c(
    "hsapiens" = "org.Hs.eg.db",
    "mmusculus" = "org.Mm.eg.db",
    "rnorvegicus" = "org.Rn.eg.db",
    "cfamiliaris" = "org.Cf.eg.db",
    "btaurus" = "org.Bt.eg.db",
    "sscrofa" = "org.Ss.eg.db",
    "ggallus" = "org.Gg.eg.db",
    "drerio" = "org.Xl.eg.db",
    "xlaevis" = "org.Dr.eg.db",
    "dmelanogaster" = "org.Dm.eg.db",
    "celegans" = "org.Ce.eg.db"
)

library(organism_mapping[args$organism], character.only = TRUE)


## ----Load omics data----------------------------------------------------------

layer <- c()

if (!is.null(args$transcriptomics)) {
    transcriptome <- read.csv(
        args$transcriptomics,
        header = TRUE,
        sep = "\t",
        dec = "."
    )
    layer <- append(layer, "transcriptome")
}

if (!is.null(args$proteomics)) {
    proteome <- read.csv(args$proteomics,
        header = TRUE,
        sep = "\t",
        dec = "."
    )
    layer <- append(layer, "proteome")
}

if (!is.null(args$metabolomics)) {
    metabolome <- read.csv(args$metabolomics,
        header = TRUE,
        sep = "\t",
        dec = "."
    )
    layer <- append(layer, "metabolome")
}

## ----rank_features------------------------------------------------------------

# create data structure
omics_data <- initOmicsDataStructure(layer)

## add transcriptome layer
if (!is.null(args$transcriptomics)) {
    omics_data$transcriptome <- rankFeatures(
        transcriptome$logFC,
        transcriptome$pValue
    )
    names(omics_data$transcriptome) <- transcriptome$Symbol
}

## add proteome layer
if (!is.null(args$proteomics)) {
    omics_data$proteome <- rankFeatures(proteome$logFC, proteome$pValue)
    names(omics_data$proteome) <- proteome$Symbol
}

## add metabolome layer
## HMDB features have to be updated to the new HMDB format
if (!is.null(args$metabolomics)) {
    omics_data$metabolome <-
        rankFeatures(metabolome$logFC, metabolome$pValue)
    names(omics_data$metabolome) <- metabolome$HMDB
    names(omics_data$metabolome) <- gsub(
        "HMDB", "HMDB00",
        names(omics_data$metabolome)
    )
}


## remove NA's and sort feature ranks
omics_data <- lapply(omics_data, function(vec) {
    sort(vec[!is.na(vec)])
})

## ----Pathway definitions------------------------------------------------------

pathways <-
    getMultiOmicsFeatures(
        dbs = unlist(strsplit(args$databases, ",", fixed = TRUE)),
        layer = layer,
        returnTranscriptome = args$transcriptome_ids,
        returnProteome = args$proteome_ids,
        returnMetabolome = args$metabolome_ids,
        organism = args$organism,
        useLocal = FALSE
    )

## ----calculate enrichment-----------------------------------------------------

enrichment_scores <-
    multiGSEA(pathways, omics_data)

## ----combine_pvalues----------------------------------------------------------

df <- extractPvalues(
    enrichmentScores = enrichment_scores,
    pathwayNames = names(pathways[[1]])
)

df$combined_pval <-
    combinePvalues(df, method = args$combine_pvalues)
df$combined_padj <-
    p.adjust(df$combined_pval, method = args$padj_method)

df <- cbind(data.frame(pathway = names(pathways[[1]])), df)

## ----Write output-------------------------------------------------------------

write.table(
    df,
    file = "results.tsv",
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
)
