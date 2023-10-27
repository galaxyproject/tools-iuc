# setup R error handling to go to stderr
options(show.error.messages = FALSE, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, FALSE)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
library("getopt")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "bamdir", "b", 1, "character",
  "gtffile", "g", 1, "character",
  "codon_coverage_info", "Y", 1, "character",
  "cds_coverage_info", "Z", 1, "character",
  "psite_info_rdata", "O", 0, "character",
  "refseq_sep", "s", 0, "character",
  "indel_threshold", "t", 0, "integer",
  "params_duplicate_filterting", "d", 0, "character",
  "params_peridiocity_filterting", "l", 0, "character",
  "params_custom_filterting", "c", 0, "character",
  "params_psite_additional", "p", 0, "character",
  "params_coverage_additional", "C", 0, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

verbose <- is.null(opt$quiet)

library("riboWaltz")

# create annotation data table
annotation_dt <- create_annotation(opt$gtffile)

sep <- opt$refseq_sep
if (opt$refseq_sep == "") {
  sep <- NULL
}
# convert alignments in BAM files into list of data tables
reads_list <- bamtolist(bamfolder = opt$bamdir, annotation = annotation_dt, refseq_sep = sep, indel_threshold = opt$indel_threshold)

library("jsonlite")
# remove duplicate reads
if (!is.null(opt$params_duplicate_filterting)) {
  json_duplicate_filterting <- fromJSON(opt$params_duplicate_filterting)
  reads_list <- duplicates_filter(
    data = reads_list,
    extremity = json_duplicate_filterting$extremity,
    keep = json_duplicate_filterting$keep
  )
}

# selection of read lengths - periodicity filtering
if (!is.null(opt$params_peridiocity_filterting)) {
  json_peridiocity_filterting <- fromJSON(opt$params_peridiocity_filterting)
  reads_list <- length_filter(
    data = reads_list,
    length_filter_mode = "periodicity",
    periodicity_threshold = json_peridiocity_filterting$periodicity_threshold
  )
}
# selection of read lengths - length range filtering
if (!is.null(opt$params_custom_filterting)) {
  json_custom_filterting <- fromJSON(opt$params_custom_filterting)
  reads_list <- length_filter(
    data = reads_list,
    length_filter_mode = "custom",
    length_range = json_custom_filterting$length_range
  )
}

# compute P-site offset
json_psite_additional <- fromJSON(opt$params_psite_additional)
psite_offset <- psite(
  reads_list,
  start = json_psite_additional$use_start,
  flanking = json_psite_additional$flanking,
  extremity = json_psite_additional$psite_extrimity,
  plot = TRUE,
  cl = json_psite_additional$cl,
  plot_format = "pdf",
  plot_dir = "plots"
)
psite_offset

reads_psite_list <- psite_info(reads_list, psite_offset)
reads_psite_list
# write a separate P-site offset info table for each sample
for (sample in names(reads_psite_list)) {
  write.table(
    reads_psite_list[[sample]],
    file = paste(sample, "psite_info.tsv",  sep = "_"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  print(paste(sample, "psite_info.tsv",  sep = "_"))
}

# write R object to a file
if (!is.null(opt$psite_info_rdata)) {
  save(reads_psite_list, annotation_dt, file = opt$psite_info_rdata)
}

json_coverage_additional <- fromJSON(opt$params_coverage_additional)
# codon coverage
codon_coverage_out <- codon_coverage(
  reads_psite_list,
  annotation_dt,
  psite = json_coverage_additional$psites_per_region,
  min_overlap = json_coverage_additional$min_overlap
)
write.table(codon_coverage_out, file = opt$codon_coverage_info, sep = "\t", row.names = FALSE, quote = FALSE)

# CDS coverage
cds_coverage_out <- cds_coverage(
  reads_psite_list,
  annotation_dt,
  start_nts = json_coverage_additional$start_nts,
  stop_nts = json_coverage_additional$stop_nts
)
write.table(cds_coverage_out, file = opt$cds_coverage_info, sep = "\t", row.names = FALSE, quote = FALSE)
