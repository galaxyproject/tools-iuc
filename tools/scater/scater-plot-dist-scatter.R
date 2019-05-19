#!/usr/bin/env Rscript

# Plot the distribution of read counts and feature counts, side by side, then a scatter plot of read counts vs feature counts below

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(ggpubr))
suppressPackageStartupMessages(require(scater))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A serialized SingleCellExperiment object file in RDS format."
  ),
  make_option(
    c("-o", "--output-plot-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path of the PDF output file to save plot to."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_plot_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Input from serialized R object

sce <- readRDS(opt$input_object_file)

#do the scatter plot of reads vs genes
total_counts <- sce$total_counts
total_features <- sce$total_features
count_feats <- cbind(total_counts, total_features)
cf_dm <- as.data.frame(count_feats)

#calculate binwidths for reads and features plots. Use 20 bins
read_bins <- max(total_counts/1e6)/20
feat_bins <- max(total_features)/20

#make the plots
plot <- ggplot(cf_dm, aes(x=total_counts/1e6, y=total_features)) + geom_point(shape=1) + geom_smooth() + xlab("Read count (millions)") +
   ylab("Feature count") + ggtitle("Scatterplot of reads vs features")
plot1 <- qplot(total_counts/1e6, geom="histogram", binwidth = read_bins, ylab="Number of cells", xlab = "Read counts (millions)", fill=I("darkseagreen3")) + ggtitle("Read counts per cell")
plot2 <- qplot(total_features, geom="histogram", binwidth = feat_bins, ylab="Number of cells", xlab = "Feature counts", fill=I("darkseagreen3")) + ggtitle("Feature counts per cell")
plot3 <- plotColData(sce, y="pct_counts_MT", x="total_features_by_counts") + ggtitle("% MT genes") + geom_point(shape=1) + theme(text = element_text(size=15)) + theme(plot.title = element_text(size=15))

final_plot <- ggarrange(plot1, plot2, plot, plot3, ncol=2, nrow=2)
ggsave(opt$output_plot_file, final_plot, device="pdf")
