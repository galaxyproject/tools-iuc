# setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr()); q("no", 1, F)
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
  "out_pdf", "P", "1", "character",
  "psite_info", "X", "1", "character",
  "codon_coverage_info", "Y", "1", "character",
  "cds_coverage_info", "Z", "1", "character",
  "refseq_sep", "s", 0, "character",
  "params_duplicate_filterting", "d", 0, "character",
  "params_peridiocity_filterting", "l", 0, "character",
  "params_custom_filterting", "c", 0, "character",
  "params_psite_additional", "p", 0, "character",
  "params_coverage_additional", "C", 0, "character",
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

library(riboWaltz)

# create annotation data table
annotation_dt <- create_annotation(opt$gtffile)

sep <- opt$refseq_sep
if (opt$refseq_sep == "") {
	sep <- NULL
}
# convert alignments in BAM files into list of data tables 
reads_list <- bamtolist(bamfolder = opt$bamdir, annotation = annotation_dt, refseq_sep=sep)

library('jsonlite')
# remove duplicate reads
if (!is.null(opt$params_duplicate_filterting)) {
  json_duplicate_filterting <- fromJSON(opt$params_duplicate_filterting)
	reads_list <- duplicates_filter(data = reads_list, extremity = json_duplicate_filterting$extremity, keep=json_duplicate_filterting$keep)
}

# selection of read lengths - periodicity filtering
if (!is.null(opt$params_peridiocity_filterting)) {
    json_peridiocity_filterting <- fromJSON(opt$params_peridiocity_filterting)
    reads_list <- length_filter(data = reads_list,
          length_filter_mode = "periodicity",
          periodicity_threshold = json_peridiocity_filterting$periodicity_threshold)
}
# selection of read lengths - length range filtering
if(!is.null(opt$params_custom_filterting)) {
    json_custom_filterting <- fromJSON(opt$params_custom_filterting)
    reads_list <- length_filter(data = reads_list,
          length_filter_mode = "custom",
          length_range=json_custom_filterting$length_range)  
}

# compute P-site offset
json_psite_additional <- fromJSON(opt$params_psite_additional)

psite_offset <- psite(reads_list, start=json_psite_additional$use_start, flanking = json_psite_additional$flanking, extremity = json_psite_additional$psite_extrimity, plot=TRUE, cl=json_psite_additional$cl, plot_format="pdf", plot_dir="plots")
psite_offset
reads_psite_list <- psite_info(reads_list, psite_offset)
reads_psite_list
json_coverage_additional <- fromJSON(opt$params_coverage_additional)
#codon_coverage_out <- codon_coverage(reads_psite_list, annotation_dt, psite = json_coverage_additional$psites_per_region, min_overlap=json_coverage_additional$min_overlap)
#cds_coverage_out <- cds_coverage(reads_psite_list, annotation_dt, start_nts=json_coverage_additional$start_nts, stop_nts=json_coverage_additional$stop_nts)


for (sample in names(reads_psite_list)) {
  write.table(reads_psite_list[[sample]], file = paste(sample, "psite_info.tsv",  sep="_"), sep = "\t", col.names = NA, quote = FALSE)
}

print(getwd())
print(paste(sample, "psite_info.tsv",  sep="_"))
#write.table(codon_coverage_out, file = opt$codon_coverage_info, sep = "\t", col.names = NA, quote = FALSE)
#write.table(cds_coverage_out, file = opt$cds_coverage_info, sep = "\t", col.names = NA, quote = FALSE)


if (!is.null(opt$params_rlength_distr) || !is.null(opt$params_rends_heat) || !is.null(opt$region_psite_plot) || !is.null(opt$params_trint_periodicity) || !is.null(opt$params_metaplots) || !is.null(opt$params_codon_usage_psite)) {
  pdf(opt$out_pdf)
}

if (!is.null(opt$params_rlength_distr)) {
  json_rlength_distr <- fromJSON(opt$params_rlength_distr)
  length_dist <- rlength_distr(reads_list, sample = names(reads_list), cl=json_rlength_distr$cl, multisamples=json_rlength_distr$multisamples, plot_style=json_rlength_distr$plot_style)
  length_dist
}

if (!is.null(opt$params_rends_heat)) {
  json_rends_heat <- fromJSON(opt$params_rends_heat)
  for (sample_name in names(reads_list)) {
	  ends_heatmap <- rends_heat(reads_list, annotation_dt, sample = sample_name, cl=json_rends_heat$cl, utr5l=json_rends_heat$utr5l, cdsl=json_rends_heat$cdsl, utr3l=json_rends_heat$utr3l)
	  ends_heatmap[["plot"]]
  }
}

if (!is.null(opt$region_psite_plot)) {
  psite_region <- region_psite(reads_psite_list, annotation_dt, sample = names(reads_list))
  psite_region[["plot"]]
}

if (!is.null(opt$params_trint_periodicity)) {
  json_trint_periodicity <- fromJSON(opt$params_trint_periodicity)
  frames_stratified <- frame_psite_length(reads_psite_list, sample = names(reads_list), cl = json_trint_periodicity$cl, region = json_trint_periodicity$region, length_range=json_trint_periodicity$length_range)
  frames_stratified[["plot"]]
  frames <- frame_psite_length(reads_psite_list, sample = names(reads_list), region = json_trint_periodicity$region, length_range=json_trint_periodicity$length_range)
  frames[["plot"]]
}

if (!is.null(opt$params_metaplots)) {
  json_metaplots <- fromJSON(opt$params_metaplots)
  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = names(reads_list), multisamples=json_metaplots$multisamples, plot_style=json_metaplots$plot_style, length_range=json_metaplots$length_range, frequency=json_metaplots$frequency, utr5l = json_metaplots$utr5l, cdsl = json_metaplots$cdsl, utr3l = json_metaplots$utr3l, plot_title = "sample.transcript.length_range")
  metaprofile
  sample_list <- list()
  for (sample_name in names(reads_list)) {
	sample_list[[sample_name]] <- c(sample_name)
  }
  metaheatmap <- metaheatmap_psite(reads_psite_list, annotation_dt, sample = sample_list, length_range=json_metaplots$length_range, utr5l = json_metaplots$utr5l, cdsl = json_metaplots$cdsl, utr3l = json_metaplots$utr3l, plot_title = "Comparison metaheatmap")
  metaheatmap[["plot"]]
}

if (!is.null(opt$params_codon_usage_psite)) {
  json_codon_usage_psite <- fromJSON(opt$params_codon_usage_psite)
  cu_barplot <- codon_usage_psite(reads_psite_list, annotation_dt, sample = names(reads_list), fastapath = json_codon_usage_psite$fastapath, fasta_genome = FALSE, frequency_normalization = json_codon_usage_psite$frequency) 
  cu_barplot[["plot"]]
}

if (!is.null(opt$params_rlength_distr) || !is.null(opt$params_rends_heat) || !is.null(opt$region_psite_plot) || !is.null(opt$params_trint_periodicity) || !is.null(opt$params_metaplots) || !is.null(opt$params_codon_usage_psite)) {
  dev.off()
}
