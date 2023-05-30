library(multiGSEA,
        quietly = TRUE,
        warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)

################################################################################
### Input Processing
################################################################################


# Collect arguments from command line
parser <- ArgumentParser(description = "multiGSEA R script")

parser$add_argument("--transcriptomics",  required = FALSE,
                    help = "Transcriptomics data")
parser$add_argument("--transcriptome_ids",  required = FALSE,
                    help = "Transcriptomics ids")
parser$add_argument("--proteomics",  required = FALSE,
                    help = "Proteomics data")
parser$add_argument("--proteome_ids",  required = FALSE,
                    help = "Proteomics ids")
parser$add_argument("--metabolomics",  required = FALSE,
                    help = "Metabolomics data")
parser$add_argument("--metabolome_ids",  required = FALSE,
                    help = "Metabolomics ids")
parser$add_argument("--organism",  required = TRUE,
                    help = "Organism")
parser$add_argument("--combine_pvalues",  required = TRUE,
                    help = "Combine p-values method")
parser$add_argument("--padj_method",  required = TRUE,
                    help = "P-adjustment method")

args <- parser$parse_args()

## --- Load library

organism_mapping <- c(
  "hsapiens"="org.Hs.eg.db",
  "mmusculus"="org.Mm.eg.db",
  "rnorvegicus"="org.Rn.eg.db",
  "cfamiliaris"="org.Cf.eg.db",
  "btaurus"="org.Bt.eg.db",
  "sscrofa"="org.Ss.eg.db",
  "ggallus"="org.Gg.eg.db",
  "drerio"="org.Xl.eg.db",
  "xlaevis"="org.Dr.eg.db",
  "dmelanogaster"="org.Dm.eg.db",
  "celegans"="org.Ce.eg.db"
)

library(organism_mapping[args$organism],character.only = TRUE)
#library(organism_mapping["hsapiens"],character.only = TRUE)

## --- Load omics data

layer <- c()

if (!is.null(args$transcriptomics)) {
  transcriptome <- read.csv(args$transcriptomics, header=TRUE,
                                   sep="\t", dec=".")
  layer <- append(layer, 'transcriptome')
}

if (!is.null(args$proteomics)) {
  proteome <- read.csv(args$proteomics, header=TRUE,
                              sep="\t", dec=".")
  layer <- append(layer, 'proteome')
}

if (!is.null(args$metabolomics)) {
  metabolome <- read.csv(args$metabolomics, header=TRUE,
                              sep="\t", dec=".")
  layer <- append(layer, 'metabolome')
}



##transcriptome <- read.csv("~/Galaxy/tools-iuc/tools/multigsea/test-data/transcriptome.tsv", header=TRUE,
##                                 sep="\t", dec=".")

##proteome <- read.csv("~/Galaxy/tools-iuc/tools/multigsea/test-data/proteome.tsv", header=TRUE,
##                                 sep="\t", dec=".")

##metabolome <- read.csv("~/Galaxy/tools-iuc/tools/multigsea/test-data/metabolome.tsv", header=TRUE,
##                                 sep="\t", dec=".")


## ----rank_features, results='hide'--------------------------------------------

# create data structure
omics_data <- initOmicsDataStructure(layer)

## add transcriptome layer
omics_data$transcriptome <- rankFeatures( transcriptome$logFC, 
                                          transcriptome$pValue)
names( omics_data$transcriptome) <- transcriptome$Symbol

## add proteome layer
omics_data$proteome <- rankFeatures(proteome$logFC, proteome$pValue)
names( omics_data$proteome) <- proteome$Symbol

## add metabolome layer
## HMDB features have to be updated to the new HMDB format
omics_data$metabolome <- rankFeatures(metabolome$logFC, metabolome$pValue)
names( omics_data$metabolome) <- metabolome$HMDB
names( omics_data$metabolome) <- gsub( "HMDB", "HMDB00", 
                                       names( omics_data$metabolome))


## ----omics_short--------------------------------------------------------------

omics_short <- lapply( names( omics_data), function( name){ 
                        head( omics_data[[name]])
                      })
names( omics_short) <- names( omics_data)

## ----calculate_enrichment, results='hide', message=FALSE, warning=FALSE-------

pathways <- getMultiOmicsFeatures( dbs = c("all"), layer = layer,
                                   returnTranscriptome = args$transcriptome_ids ,
                                   returnProteome = args$proteome_ids,
                                   returnMetabolome = args$metabolome_ids,
                                   organism = args$organism,
                                   useLocal = FALSE)

## ----pathways_short-----------------------------------------------------------

pathways_short <- lapply( names( pathways), function( name){
                          head( pathways[[name]], 2)
                        })
names( pathways_short) <- names( pathways)
pathways_short

## ----run_enrichment, results='hide', message=FALSE, warning=FALSE-------------

# use the multiGSEA function to calculate the enrichment scores

enrichment_scores <- multiGSEA( pathways, omics_data)


## ----combine_pvalues----------------------------------------------------------

df <- extractPvalues( enrichmentScores = enrichment_scores,
                      pathwayNames = names( pathways[[1]]))

df$combined_pval <- combinePvalues( df, method = args$combine_pvalues)
df$combined_padj <- p.adjust( df$combined_pval, method = args$padj_method)

df <- cbind( data.frame( pathway = names( pathways[[1]])), df)

## ----session, echo=FALSE------------------------------------------------------
sessionInfo()

