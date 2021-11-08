suppressPackageStartupMessages(library(optparse))

### Parsing command line -------------------------------------------------------
option_list <- list(
    make_option(c("--input"), type = "character",
                default = NULL, help = "Count matrix in tsv format"),
    make_option(c("--model"), type = "character",
                default = NULL, help = "RAVmodel to apply.
                Currently 'C2' and 'PLIERpriors' are available"),
    make_option(c("--method"), type = "character",
                default = formals(GenomicSuperSignature::validate)$method),
    make_option(c("--maxFrom"), type = "character",
                default = formals(GenomicSuperSignature::validate)$maxFrom),
    make_option(c("--level"), type = "character",
                default = formals(GenomicSuperSignature::validate)$level),
    make_option(c("--scale"), type = "character",
                default = formals(GenomicSuperSignature::validate)$scale),
    make_option(c("--outDir"), type = "character",
                default = NULL, help = "Output file name"),
    make_option(c("--numOut"), type = "integer",
                default = 3, help = "The number of top validated RAVs to check")
)

opt <- parse_args(OptionParser(option_list = option_list),
                  args = commandArgs(trailingOnly = TRUE))

if (is.null(opt$input)) stop("Need --input.")
if (is.null(opt$model)) stop("Need --model.")
if (is.null(opt$outDir)) stop("Need --outDir.")

inputName <- tools::file_path_sans_ext(opt$input)

suppressPackageStartupMessages(library(GenomicSuperSignature))
dat <- as.matrix(read.table(file = opt$input, header = TRUE, sep = "\t"))
# RAVmodel <- getModel(opt$model)
RAVmodel <! readRDS(opt$model)


### validate -------------------------------------------------------------------
val_all <- validate(dat, RAVmodel)
output_fname <- paste0(inputName, "_validate.csv")

write.csv(val_all,
          file = file.path(opt$outDir, output_fname),
          row.names = TRUE)

### Save graphs in pdf ---------------------------------------------------------
pdf("gss_outputs.pdf")

# heatmap table of validation result
heatmapTable(val_all, num.out = opt$numOut, swCutoff = 0)

# interactive plot for validation result
plotValidate(val_all, interactive = FALSE)

# MeSH term wordcloud
validated_ind <- validatedSignatures(val_all, num.out = opt$numOut,
                                     swCutoff = 0, indexOnly = TRUE)
for (i in seq_len(opt$numOut)) {
    set.seed(1)
    drawWordcloud(RAVmodel, validated_ind[i])
}
dev.off()

### Save tables in csv ---------------------------------------------------------
# GSEA
for (i in seq_len(opt$numOut)) {
    RAVnum <- validated_ind[i]
    res <- gsea(RAVmodel)[[RAVnum]]

    output_fname <- paste0(inputName, "_genesets_RAV", RAVnum, ".csv")
    write.csv(res,
              file = file.path(opt$outDir, output_fname),
              row.names = TRUE)
}

# Related prior studies
for (i in seq_len(opt$numOut)) {
    RAVnum <- validated_ind[i]
    findStudiesInCluster(RAVmodel, RAVnum)

    output_fname <- paste0(inputName, "_literatures_RAV", RAVnum, ".csv")
    write.csv(res,
              file = file.path(opt$outDir, output_fname),
              row.names = TRUE)
}


