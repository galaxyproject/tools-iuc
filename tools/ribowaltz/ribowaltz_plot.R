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
  "input_rdata", "i", 1, "character",
  "params_rlength_distr", "r", 0, "character",
  "params_rends_heat", "e", 0, "character",
  "region_psite_plot", "R", 0, "logical",
  "params_trint_periodicity", "t", 0, "character",
  "params_metaplots", "m", 0, "character",
  "params_codon_usage_psite", "u", 0, "character"
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
library("jsonlite")

load(opt$input_rdata)

if (!is.null(opt$params_rlength_distr)) {
  pdf("read_lengths.pdf")
  json_rlength_distr <- fromJSON(opt$params_rlength_distr)
  length_dist <- rlength_distr(
    reads_psite_list,
    sample = names(reads_psite_list),
    cl = json_rlength_distr$cl,
    multisamples = json_rlength_distr$multisamples,
    plot_style = json_rlength_distr$plot_style
  )
  print(length_dist)
  dev.off()
}

if (!is.null(opt$params_rends_heat)) {
  pdf("read_ends_heatmap.pdf", height = 5 * length(reads_psite_list), width = 15)
  json_rends_heat <- fromJSON(opt$params_rends_heat)
  for (sample_name in names(reads_psite_list)) {
    ends_heatmap <- rends_heat(
      reads_psite_list,
      annotation_dt,
      sample = sample_name,
      cl = json_rends_heat$cl,
      utr5l = json_rends_heat$utr5l,
      cdsl = json_rends_heat$cdsl,
      utr3l = json_rends_heat$utr3l
    )
    print(ends_heatmap[["plot"]])
  }
  dev.off()
}

if (!is.null(opt$region_psite_plot)) {
  pdf("psites_per_region.pdf", height = 12, width = 7 * length(reads_psite_list))
  psite_region <- region_psite(reads_psite_list, annotation_dt, sample = names(reads_psite_list))
  print(psite_region[["plot"]])
  dev.off()
}

if (!is.null(opt$params_trint_periodicity)) {
  pdf("trinucleotide_periodicity.pdf", height = 6 * length(reads_psite_list), width = 10)
  json_trint_periodicity <- fromJSON(opt$params_trint_periodicity)
  frames_stratified <- frame_psite_length(
    reads_psite_list,
    sample = names(reads_psite_list),
    cl = json_trint_periodicity$cl,
    region = json_trint_periodicity$region,
    length_range = json_trint_periodicity$length_range
  )
  frames_stratified[["plot"]]
  frames <- frame_psite_length(
    reads_psite_list,
    sample = names(reads_psite_list),
    region = json_trint_periodicity$region,
    length_range = json_trint_periodicity$length_range
  )
  print(frames[["plot"]])
  dev.off()
}

if (!is.null(opt$params_metaplots)) {
  pdf("metaplots.pdf", height = 5 * length(reads_psite_list), width = 24)
  json_metaplots <- fromJSON(opt$params_metaplots)
  metaprofile <- metaprofile_psite(
    reads_psite_list,
    annotation_dt,
    sample = names(reads_psite_list),
    multisamples = json_metaplots$multisamples,
    plot_style = json_metaplots$plot_style,
    length_range = json_metaplots$length_range,
    frequency = json_metaplots$frequency,
    utr5l = json_metaplots$utr5l,
    cdsl = json_metaplots$cdsl,
    utr3l = json_metaplots$utr3l,
    plot_title = "sample.transcript.length_range"
  )
  print(metaprofile)
  sample_list <- list()
  for (sample_name in names(reads_psite_list)) {
    sample_list[[sample_name]] <- c(sample_name)
  }
  metaheatmap <- metaheatmap_psite(
    reads_psite_list,
    annotation_dt,
    sample = sample_list,
    length_range = json_metaplots$length_range,
    utr5l = json_metaplots$utr5l,
    cdsl = json_metaplots$cdsl,
    utr3l = json_metaplots$utr3l,
    plot_title = "Comparison metaheatmap"
  )
  print(metaheatmap[["plot"]])
  dev.off()
}

if (!is.null(opt$params_codon_usage_psite)) {
  pdf("codon_usage.pdf", height = 6, width = 16)
  json_codon_usage_psite <- fromJSON(opt$params_codon_usage_psite)
  for (sample_name in names(reads_psite_list)) {
    cu_barplot <- codon_usage_psite(
    reads_psite_list,
    annotation_dt,
    sample = sample_name,
    fastapath = json_codon_usage_psite$fastapath,
    fasta_genome = FALSE,
    frequency_normalization = json_codon_usage_psite$frequency
  )
  print(cu_barplot[["plot"]])
  }
  dev.off()
}
