library("getopt")
suppressPackageStartupMessages(library("rtracklayer"))
library(GenomicRanges)
library("BREW3R.r")

options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)
# - Column 1: the long flag name. A multi-character string.
# - Column 2: short flag alias of Column 1. A single-character string.
# - Column 3: Argument mask of the flag. An integer.
# Possible values: 0=no argument, 1=required argument, 2=optional argument.
# - Column 4: Data type to which the flag's argument shall be cast using
# storage.mode(). A multi-character string. This only considered for same-row
# Column 3 values of 1,2. Possible values: logical, integer, double, complex,
# character. If numeric is encountered then it will be converted to double.
# - Column 5 (optional): A brief description of the purpose of the option.
spec <- matrix(c(
    "help", "h", 0, "logical", "display help",
    "gtf_to_extend", "i", 1, "character", "input gtf file to be extended on 3'",
    "gtf_to_overlap", "g", 1, "character",
    "input gtf file that will be used to extend",
    "output", "o", 1, "character", "output extended gtf",
    "sup_output", "s", 1, "character",
    "supplementary output file with resolution of overlaps",
    "no_add", "n", 0, "logical", "do not add new exons",
    "exclude_pattern", "e", 1, "character", "do not extend genes with names matching this pattern",
    "filter_unstranded", "f", 0, "logical",
    "remove unstranded intervals from gtf_to_overlap which overlap intervals from gtf_to_extend of both strands",
    "quiet", "q", 0, "logical", "decrease verbosity",
    "verbose", "v", 0, "logical", "increase verbosity"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

# Check all required arguments
if (is.null(opt$gtf_to_extend)) {
    stop("--gtf_to_extend is required")
}
if (is.null(opt$gtf_to_overlap)) {
    stop("--gtf_to_overlap is required")
}
if (is.null(opt$output)) {
    stop("--output is required")
}

# Check incompatible arguments
if (!is.null(opt$quiet) && !is.null(opt$verbose)) {
    stop("quiet and verbose are mutually exclusive options")
}

# Adjust verbosity
if (!is.null(opt$quiet)) {
    options(rlib_message_verbosity = "quiet")
}

if (!is.null(opt$verbose)) {
    options(BREW3R.r.verbose = "progression")
}

# Load gtfs as GenomicRanges
input_gr_to_extend <- rtracklayer::import(opt$gtf_to_extend, format = "gtf")
input_gr_template <- rtracklayer::import(opt$gtf_to_overlap, format = "gtf")

# Save CDS info
input_gr_CDS <- subset(input_gr_to_extend, type == "CDS")

# Filter the template if needed
if (!is.null(opt$filter_unstranded)) {
    # Find intervals without strand information in template
    unstranded.intervals <- which(strand(input_gr_template) == "*")
    if (length(unstranded.intervals) > 0) {
        # Check if they overlap genes from input with different strands
        # First compute the overlap
        ov <- suppressWarnings(
            as.data.frame(findOverlaps(
                input_gr_template[unstranded.intervals],
                input_gr_to_extend
            ))
        )
        # Add the strand information
        ov$strand <- as.factor(strand(input_gr_to_extend))[ov$subjectHits]
        # Simplify the dataframe to get only the strand info
        ov.simple <- unique(ov[, c("queryHits", "strand")])
        # If the queryHits is duplicated it means there are different strands
        multi.strand.query <- ov.simple$queryHits[duplicated(ov.simple$queryHits)]
        to.remove <- unstranded.intervals[multi.strand.query]
        # Remove these potentially error-prone intervals from the template
        if (length(to.remove) > 0) {
            input_gr_template <- input_gr_template[-to.remove]
        }
    }
}

if (is.null(input_gr_to_extend$exon_id)) {
    is.exon <- which(input_gr_to_extend$type == "exon")
    input_gr_to_extend$exon_id <- NA
    input_gr_to_extend$exon_id[is.exon] <- paste0(
        "EXON",
        sprintf(
            "%010d",
            1:length(is.exon)
        )
    )
}

# Run BREW3R.r main function
if (length(input_gr_template) > 0) {
    new_gr_exons <- extend_granges(
        input_gr_to_extend = input_gr_to_extend,
        input_gr_to_overlap = input_gr_template,
        add_new_exons = is.null(opt$no_add),
        overlap_resolution_fn = opt$sup_output
    )
} else {
    new_gr_exons <- subset(input_gr_to_extend, type == "exon")
}
# Prevent extension using pattern
if (!is.null(opt$exclude_pattern)) {
    input_gr_pattern <- subset(
        input_gr_to_extend,
        type == "exon" & grepl(opt$exclude_pattern, gene_name)
    )
    new_gr_no_pattern <- subset(
        new_gr_exons,
        !grepl(opt$exclude_pattern, gene_name)
    )
    new_gr_exons <- c(new_gr_no_pattern, input_gr_pattern)
}

# Recompose with CDS
new_gr <- c(new_gr_exons, input_gr_CDS)

# Export
rtracklayer::export.gff(sort(new_gr, ignore.strand = TRUE), opt$output)
