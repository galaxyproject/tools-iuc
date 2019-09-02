#!/usr/bin/env Rscript

# Plot the distribution of read counts and feature counts, side by side, then a scatter plot of read counts vs feature counts below

# Load optparse we need to check inputs

library(optparse)
library(workflowscriptscommon)
library(LoomExperiment)
library(scater)
library(ggpubr)
library(scales)

# parse options

option_list = list(
  make_option(
    c("-i", "--input-loom"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A SingleCellExperiment object file in Loom format."
  ),
  make_option(
    c("-o", "--output-plot-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path of the PDF output file to save plot to."
  ),
  make_option(
    c("-l", "--log-scale"),
    action="store", 
    default="FALSE",
    type = 'character',
    help = "Plot log scale (reccommended for large datasets)."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_loom', 'output_plot_file', 'log_scale'))

# Check parameter values

if ( ! file.exists(opt$input_loom)){
  stop((paste('File', opt$input_loom, 'does not exist')))
}

# Input from Loom format

scle <- import(opt$input_loom, format='loom', type='SingleCellLoomExperiment')

#do the scatter plot of reads vs genes
total_counts <- scle$total_counts
total_features <- scle$total_features_by_counts
count_feats <- cbind(total_counts, total_features)
cf_dm <- as.data.frame(count_feats)

# Calculate binwidths for reads and features plots. Use 20 bins
read_bins <- max(total_counts / 1e6) / 20
feat_bins <- max(total_features) / 20


#make the plots
plot <- ggplot(cf_dm, aes(x=total_counts / 1e6, y=total_features)) + geom_point(shape=1) + geom_smooth() + xlab("Read count (millions)") +
   ylab("Feature count") + ggtitle("Scatterplot of reads vs features")
plot1 <- qplot(total_counts / 1e6, geom="histogram", binwidth = read_bins, ylab="Number of cells", xlab = "Read counts (millions)", fill=I("darkseagreen3")) + ggtitle("Read counts per cell")
plot2 <- qplot(total_features, geom="histogram", binwidth = feat_bins, ylab="Number of cells", xlab = "Feature counts", fill=I("darkseagreen3")) + ggtitle("Feature counts per cell")
plot3 <- plotColData(scle, y="pct_counts_MT", x="total_features_by_counts") + ggtitle("% MT genes") + geom_point(shape=1) + theme(text = element_text(size=15)) + theme(plot.title = element_text(size=15))


if (identical(opt$log_scale, 'FALSE'))
{
  final_plot <- ggarrange(plot1, plot2, plot, plot3, ncol=2, nrow=2)
  ggsave(opt$output_plot_file, final_plot, device="pdf")
} else 
  {
    plot_log_both <- plot + scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10')
    plot1_log <- plot1 + scale_y_continuous(trans = 'log10')
    plot2_log <- plot2 + scale_y_continuous(trans = 'log10')
    plot3_log <- plot3 + scale_y_log10(labels=number)
    final_plot_log <- ggarrange(plot1_log, plot2_log, plot_log_both, plot3_log, ncol=2, nrow=2)
    ggsave(opt$output_plot_file, final_plot_log, device="pdf")
}