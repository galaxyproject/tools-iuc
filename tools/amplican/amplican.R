## Setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr()); q("no", 1, F)
})
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("amplican")
    library("getopt")
})

# Collect arguments from command line
args <- commandArgs(trailingOnly = TRUE)
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "config", "c", 1, "character",
  "fastq_folder", "i", 1, "character",
  "results_folder", "o", 1, "character",
  "knit_reports", "k", 1, "logical",
  "write_alignments_format", "f", 1, "character",
  "average_quality", "a", 1, "integer",
  "min_quality", "m", 1, "integer",
  "use_parallel", "u", 1, "logical",
  "match_scoring", "M", 1, "integer",
  "mismatch_scoring", "s", 1, "integer",
  "base_only", "b", 1, "logical",
  "type", "t", 1, "character",
  "gap_opening", "O", 1, "character",
  "gap_extension", "e", 1, "character",
  "fastqfiles", "T", 1, "character",
  "primer_mismatch", "P", 1, "integer",
  "donor_mismatch", "A", 1, "integer",
  "primer_dimer", "B", 1, "integer",
  "event_filter", "F", 1, "logical",
  "cut_buffer", "d", 1, "integer",
  "promiscuous_consensus", "C", 1, "logical",
  "min_freq", "Q", 1, "numerical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# path to example config file
config <- system.file("extdata", "config.csv", package = "amplican")
# path to example fastq files
fastq_folder <- system.file("extdata", package = "amplican")
# output folder
results_folder <- tempdir()

#full analysis, not knitting files automatically
amplicanPipeline(opt$config, 
  opt$fastq_folder,
  opt$results_folder,
  knit_reports = opt$knit_reports,
  write_alignments_format = c(opt$write_alignments_format),
  average_quality = opt$average_quality,
  min_quality = opt$min_quality,
  use_parallel = opt$use_parallel,
  scoring_matrix = Biostrings::nucleotideSubstitutionMatrix(match = opt$match_scoring,
    mismatch = opt$mismatch_scoring, 
    baseOnly = opt$base_only, 
    type = opt$type),
  gap_opening = opt$gap_opening,
  gap_extension = opt$gap_extension, 
  fastqfiles = opt$fastqfiles, 
  primer_mismatch = opt$primer_mismatch,
  donor_mismatch = opt$donor_mismatch, 
  PRIMER_DIMER = opt$primer_dimer, 
  event_filter = opt$event_filter,
  cut_buffer = opt$cut_buffer, 
  promiscuous_consensus = opt$promiscuous_consensus,
  normalize = c("guideRNA", "Group"), 
  min_freq = opt$min_freq)