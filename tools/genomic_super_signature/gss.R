suppressPackageStartupMessages(library(optparse))

### Parsing command line -------------------------------------------------------
option_list <- list(
    make_option(c("--input"),
        type = "character",
        default = NULL, help = "Count matrix in tsv format"
    ),
    make_option(c("--model"),
        type = "character",
        default = NULL, help = "RAVmodel to apply.
                Currently 'C2' and 'PLIERpriors' are available"
    ),
    make_option(c("--method"),
        type = "character",
        default = formals(GenomicSuperSignature::validate)$method
    ),
    make_option(c("--maxFrom"),
        type = "character",
        default = formals(GenomicSuperSignature::validate)$maxFrom
    ),
    make_option(c("--level"),
        type = "character",
        default = formals(GenomicSuperSignature::validate)$level
    ),
    make_option(c("--scale"),
        type = "character",
        default = formals(GenomicSuperSignature::validate)$scale
    ),
    make_option(c("--outDir"),
        type = "character",
        default = NULL, help = "Output file name"
    ),
    make_option(c("--validate"),
        type = "character",
        default = NULL, help = "Path to save validate.csv"
    ),
    make_option(c("--html"),
        type = "character",
        default = NULL, help = "Path to save HTML report"
    ),
    make_option(c("--numOut"),
        type = "integer",
        default = 3, help = "The number of top validated RAVs to check"
    ),
    make_option(c("--toolDir"),
        type = "character",
        default = ".", help = "Directory containing the tool scripts (e.g. gss.Rmd"
    )
)

opt <- parse_args(OptionParser(option_list = option_list),
    args = commandArgs(trailingOnly = TRUE)
)
input <- opt$input
model <- opt$model
out_dir <- opt$outDir
num_out <- opt$numOut

if (is.null(input)) stop("Need --input.")
if (is.null(model)) stop("Need --model.")
if (is.null(out_dir)) stop("Need --outDir.")

input_name <- basename(tools::file_path_sans_ext(input))
out_dir <- normalizePath(out_dir)

suppressPackageStartupMessages(library(GenomicSuperSignature))
dat <- as.matrix(read.table(
    file = input, header = TRUE, sep = "\t",
    row.names = 1
))
if (model %in% c("C2", "PLIERpriors")) {
    rav_model <- getModel(model)
} else {
    rav_model <- readRDS(model)
}



### validate -------------------------------------------------------------------
val_all <- validate(dat, rav_model)
validated_ind <- validatedSignatures(val_all,
    num.out = num_out,
    swCutoff = 0, indexOnly = TRUE
)
n <- min(num_out, length(validated_ind), na.rm = TRUE)

### Save tables in csv ---------------------------------------------------------
# Validation
if (is.null(opt$validate)) {
    output_fname <- file.path(out_dir, paste0(input_name, "_validate.csv"))
} else {
    output_fname <- opt$validate
}
write.csv(val_all,
    file = output_fname,
    row.names = TRUE
)

# GSEA
for (i in seq_len(n)) {
    rav_num <- validated_ind[i]
    rav_name <- paste0("RAV", rav_num)
    res <- gsea(rav_model)[[rav_name]]

    output_fname <- paste0(input_name, "_genesets_RAV", rav_num, ".csv")
    write.csv(res,
        file = file.path(out_dir, output_fname),
        row.names = TRUE
    )
}

# Related prior studies
for (i in seq_len(n)) {
    rav_num <- validated_ind[i]
    res <- findStudiesInCluster(rav_model, rav_num)

    output_fname <- paste0(input_name, "_literatures_RAV", rav_num, ".csv")
    write.csv(res,
        file = file.path(out_dir, output_fname),
        row.names = TRUE
    )
}

### Create a report ------------------------------------------------------------
if (is.null(opt$html)) {
    output_fname <- file.path(out_dir, paste0(
        "GSS-", input_name, "-",
        format(Sys.Date(), format = "%Y%m%d"), ".html"
    ))
} else {
    output_fname <- opt$html
}
rmarkdown::render(
    file.path(opt$toolDir, "gss.Rmd"),
    params = list(
        val_all = val_all,
        dat = dat,
        RAVmodel = rav_model,
        inputName = input_name,
        numOut = num_out
    ),
    output_file = output_fname,
    intermediates_dir = ".",
    knit_root_dir = "."
)
