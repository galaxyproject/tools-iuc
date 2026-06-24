library(SCENIC)
library(SCopeLoomR)
library(optparse)

# Set up command line options
option_list <- list(
  make_option("-i", "--input-file", dest="input_file", type="character", 
              help="Input expression matrix file (in tabular format)"),
  make_option("-m", "--motifs-file", dest="motifs_file", type="character",
              help="File containing the motifs (in tabular format)"),
  make_option("-o", "--output-dir", dest="output_dir", type="character", 
              help="Output directory"),
  make_option("-c", "--num-cores", dest="num_cores", type="integer", default=1,
              help="Number of cores to use for parallel processing (default is 1)")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# Load input expression matrix
expr_matrix <- read.table(args$input_file, header=TRUE, row.names=1, sep="\t")

# Load motifs file
motifs <- read.table(args$motifs_file, header=TRUE, sep="\t")

# Convert motifs to REGULONS object
motif_map <- create_motif_map(motifs)
regulons <- infer_regulons(expr_matrix, motif_map, num_cores=args$num_cores)

# Save regulons object to loom file
regulons_file <- file.path(args$output_dir, "regulons.loom")
write_scopeloom(regulons, file=regulons_file)

# Save correlation matrix to tsv file
cor_matrix <- compute_correlation_matrix(expr_matrix, regulons)
cor_matrix_file <- file.path(args$output_dir, "cor_matrix.tsv")
write.table(cor_matrix, file=cor_matrix_file, quote=FALSE, sep="\t", col.names=NA)

# Save enrichment scores to tsv file
enrichment_scores <- compute_enrichment_scores(expr_matrix, regulons)
enrichment_scores_file <- file.path(args$output_dir, "enrichment_scores.tsv")
write.table(enrichment_scores, file=enrichment_scores_file, quote=FALSE, sep="\t", col.names=NA)
