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
    make_option(c("--validate"), type = "character",
                default = NULL, help = "Path to save validate.csv"),
    make_option(c("--html"), type = "character",
                default = NULL, help = "Path to save HTML report"),
    make_option(c("--numOut"), type = "integer",
                default = 3, help = "The number of top validated RAVs to check"),
    make_option(c("--toolDir"), type = "character",
                default = '.', help = "Directory containing the tool scripts (e.g. gss.Rmd")
)

opt <- parse_args(OptionParser(option_list = option_list),
                  args = commandArgs(trailingOnly = TRUE))
input <- opt$input
model <- opt$model
outDir <- opt$outDir
numOut <- opt$numOut

if (is.null(input)) stop("Need --input.")
if (is.null(model)) stop("Need --model.")
if (is.null(outDir)) stop("Need --outDir.")

inputName <- basename(tools::file_path_sans_ext(input))
outDir <- normalizePath(outDir)

suppressPackageStartupMessages(library(GenomicSuperSignature))
dat <- as.matrix(read.table(file = input, header = TRUE, sep = "\t",
                            row.names = 1))
if (model %in% c("C2", "PLIERpriors")) {
    RAVmodel <- getModel(model)
} else {
    RAVmodel <- readRDS(model)
}



### validate -------------------------------------------------------------------
val_all <- validate(dat, RAVmodel)
validated_ind <- validatedSignatures(val_all, num.out = numOut,
                                     swCutoff = 0, indexOnly = TRUE)
n <- min(numOut, length(validated_ind), na.rm = TRUE)

### Save interactive plot in html ----------------------------------------------
# # interactive plot for validation result
# plot_val <- plotValidate(val_all, interactive = TRUE)
# output_fname <- paste0(inputName, "_validate_plot.html")
# htmltools::save_html(plot_val, file = file.path(outDir, output_fname))

### Save tables in csv ---------------------------------------------------------
# Validation
if (is.null(opt$validate)) {
    output_fname <- file.path(outDir, paste0(inputName, "_validate.csv"))
} else {
    output_fname <- opt$validate
}
write.csv(val_all,
          file = output_fname,
          row.names = TRUE)

# GSEA
for (i in seq_len(n)) {
    RAVnum <- validated_ind[i]
    RAVname <- paste0("RAV", RAVnum)
    res <- gsea(RAVmodel)[[RAVname]]

    output_fname <- paste0(inputName, "_genesets_RAV", RAVnum, ".csv")
    write.csv(res,
              file = file.path(outDir, output_fname),
              row.names = TRUE)
}

# Related prior studies
for (i in seq_len(n)) {
    RAVnum <- validated_ind[i]
    res <- findStudiesInCluster(RAVmodel, RAVnum)

    output_fname <- paste0(inputName, "_literatures_RAV", RAVnum, ".csv")
    write.csv(res,
              file = file.path(outDir, output_fname),
              row.names = TRUE)
}

### Create a report ------------------------------------------------------------
if (is.null(opt$html)) {
    output_fname <- file.path(outDir, paste0("GSS-", inputName, "-",
                              format(Sys.Date(), format="%Y%m%d"), ".html"))
} else {
    output_fname <- opt$html

}
rmarkdown::render(
    file.path(opt$toolDir, "gss.Rmd"), params = list(
        val_all = val_all,
        dat = dat,
        RAVmodel = RAVmodel,
        inputName = inputName,
        numOut = numOut
    ),
    output_file = output_fname
)
