## How to execute this tool

# A command-line interface to lineagespot for use with Galaxy
#
# The following arguments are required:
#
#   'in_vcf'  a character vector of paths to VCF files object from Galaxy lineagespot/test-data 
#   'in_gff3' a character vector of path to GFF3 file containing SARS-CoV-2 gene coordinates object from Galaxy lineagespot/test-data 
#   'in_ref'  a character vector of path to a folder containing lineage reports object from Galaxy lineagespot/test-data 
#   'in_voc'  a character vector containing the names of the lineages of interest 
#   'in_threshold' a parameter indicating the AF threshold for identifying variants per sample
#
#Rscript ${__tool_directory__}/lineagespot_verbose.R --in_vcf ${__tool_directory__}/test-data/extdata/vcf-files --in_gff3 ${__tool_directory__}/test-data/extdata/NC_045512.2_annot.gff3 --in_ref ${__tool_directory__}/test-data/extdata/ref --in_voc "B.1.617.2, B.1.1.7, B.1.351, P.1" --in_threshold 0.8
# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)}) 

# Import required libraries

library.path <- .libPaths()

suppressPackageStartupMessages({
  library("getopt", lib.loc = library.path)
  library("data.table", lib.loc = library.path)
  library("lineagespot", lib.loc = library.path)
}) 


options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification <- matrix(c(
  'in_vcf', 'vcf',     1, 'character',
  'in_gff3', 'gff3',   1, 'character',
  'in_ref', 'ref',     1, 'character',
  'in_voc', 'voc',     2, "character",
  'in_threshold', "thr", 2, 'double'
  # voc = c("B.1.617.2", "B.1.1.7", "B.1.351", "P.1"),
  # AF_threshold = 0.8
), byrow=TRUE, ncol =4);


# Parse options
options <- getopt(option_specification);

# Print arguments
# cat("\n vcf_file: ",options$in_vcf)
# cat("\n ref_file: ",options$in_ref)
# cat("\n gff3_file: ",options$in_gff3)

if (!is.null(options$in_voc)&is.character(options$in_voc)){
    options$in_voc = unlist(strsplit(options$in_voc,split = ','))
   
}

result <- lineagespot(vcf_folder = options$in_vcf,
                     ref_folder = options$in_ref,
                     gff3_path = options$in_gff3,
                     voc =  options$in_voc,
                     AF_threshold = options$in_threshold)


# Write output to new file which will be recognized by Galaxy
fwrite(result$variants.table, sep = '\t',file="variants_table.txt",row.names = FALSE)
fwrite(result$lineage.hits, sep = '\t', file="lineage_hits.txt",row.names = FALSE)
fwrite(result$lineage.report, sep = '\t', file="lineage_report.txt",row.names = FALSE)

cat("\n Process has been completed !\n")
