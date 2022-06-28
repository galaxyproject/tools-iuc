#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ALDEx2"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("qgraph"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--aldex_test"), action = "store", dest = "aldex_test", default = NULL, help = "Indicates which analysis to perform"),
    make_option(c("--analysis_type"), action = "store", dest = "analysis_type", help = "Indicates which analysis to perform"),
    make_option(c("--cutoff_effect"), action = "store", dest = "cutoff_effect", type = "integer", default = NULL, help = "Effect size cutoff for plotting"),
    make_option(c("--cutoff_pval"), action = "store", dest = "cutoff_pval", type = "double", default = NULL, help = "Benjamini-Hochberg fdr cutoff"),
    make_option(c("--denom"), action = "store", dest = "denom", help = "Indicates which features to retain as the denominator for the Geometric Mean calculation"),
    make_option(c("--effect"), action = "store", dest = "effect", default = "false", help = "Calculate abundances and effect sizes"),
    make_option(c("--feature_name"), action = "store", dest = "feature_name", default = NULL, help = "Name of the feature from the input data to be plotted"),
    make_option(c("--group_names"), action = "store", dest = "group_names", help = "Group names vector"),
    make_option(c("--group_nums"), action = "store", dest = "group_nums", default = NULL, help = "Group number for continuous numeric vector"),
    make_option(c("--hist_plot"), action = "store", dest = "hist_plot", default = "false", help = "Indicates whether to plot a histogram of p-values for the first Dirichlet Monte Carlo instance"),
    make_option(c("--include_sample_summary"), action = "store", dest = "include_sample_summary", default = "false", help = "Include median clr values for each sample"),
    make_option(c("--iterate"), action = "store", dest = "iterate", default = "false", help = "Indicates whether to iteratively perform a test"),
    make_option(c("--num_cols"), action = "store", dest = "num_cols", help = "Number of columns in group vector"),
    make_option(c("--num_cols_in_groups"), action = "store", dest = "num_cols_in_groups", default = NULL, help = "Number of columns in each group dewfining the continuous numeric vector"),
    make_option(c("--num_mc_samples"), action = "store", dest = "num_mc_samples", type = "integer", help = "Number of Monte Carlo samples"),
    make_option(c("--output"), action = "store", dest = "output", help = "output file"),
    make_option(c("--paired_test"), action = "store", dest = "paired_test", default = "false", help = "Indicates whether to do paired-sample tests"),
    make_option(c("--plot_test"), action = "store", dest = "plot_test", default = NULL, help = "The method of calculating significance"),
    make_option(c("--plot_type"), action = "store", dest = "plot_type", default = NULL, help = "The type of plot to be produced"),
    make_option(c("--reads"), action = "store", dest = "reads", help = "Input reads table"),
    make_option(c("--xlab"), action = "store", dest = "xlab", default = NULL, help = "x lable for the plot"),
    make_option(c("--ylab"), action = "store", dest = "ylab", default = NULL, help = "y lable for the plot")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

get_boolean_value <- function(val) {
    if (val == "true") {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Read the input reads file into a data frame.
reads_df <- read.table(file = opt$reads, header = TRUE, sep = "\t", row.names = 1, dec = ".", as.is = FALSE)

# Split the group_names and num_cols into lists of strings.
group_names_str <- as.character(opt$group_names)
group_names_list <- strsplit(group_names_str, ",")[[1]]
num_cols_str <- as.character(opt$num_cols)
num_cols_list <- strsplit(num_cols_str, ",")[[1]]
# Construct conditions vector.
conditions_vector <- c()
for (i in seq_along(num_cols_list)) {
    num_cols <- as.integer(num_cols_list[i])
    group_name <- group_names_list[i]
    for (j in 1:num_cols) {
        conditions_vector <- cbind(conditions_vector, group_name)
    }
}
# The conditions_vector is now a matrix,
# so coerce it back to a vector.
conditions_vector <- as.vector(conditions_vector)

# Convert boolean values to boolean.
effect <- get_boolean_value(opt$effect)
include_sample_summary <- get_boolean_value(opt$include_sample_summary)
iterate <- get_boolean_value(opt$iterate)

if (opt$analysis_type == "aldex") {
    aldex_obj <- aldex(reads = reads_df,
                       conditions_vector,
                       mc.samples = opt$num_mc_samples,
                       test = opt$aldex_test,
                       effect = effect,
                       include.sample.summary = include_sample_summary,
                       denom = opt$denom,
                       iterate = iterate)
} else {
    # Generate Monte Carlo samples of the Dirichlet distribution for each sample. Convert each
    # instance using a log-ratio transform. This is the input for all further analyses.
    aldex_clr_obj <- aldex.clr(reads_df, conditions_vector, mc.samples = opt$num_mc_samples, denom = opt$denom)

    if (opt$analysis_type == "aldex_corr") {
        if (!is.null(opt$cont_var)) {
            # Read the input cont_var vector.
            cont_var <- as.numeric(read.table(file = opt$cont_var, header = FALSE, sep = "\t"))
        }

        # Split the group_names and num_cols into lists of strings.
        group_nums_str <- as.character(opt$group_nums)
        group_nums_list <- strsplit(group_nums_str, ",")[[1]]
        num_cols_in_groups_str <- as.character(opt$num_cols_in_groups)
        num_cols_in_groups_list <- strsplit(num_cols_in_groups_str, ",")[[1]]
        # Construct continuous numeric vector.
        cont_var_vector <- c()
        for (i in seq_along(num_cols_in_groups_list)) {
            num_cols_in_group <- as.integer(num_cols_in_groups_list[i])
            group_num <- group_nums_list[i]
            for (j in 1:num_cols_in_group) {
                cont_var_vector <- cbind(cont_var_vector, group_num)
            }
        }
        # The cont_var_vector is now a matrix,
        # so coerce it back to a vector.
        cont_var_vector <- as.numeric(as.vector(cont_var_vector))

        aldex_obj <- aldex.corr(aldex_clr_obj, cont.var = cont_var_vector)
    } else if (opt$analysis_type == "aldex_effect") {
        aldex_obj <- aldex.effect(aldex_clr_obj, include_sample_summary)
    } else if (opt$analysis_type == "aldex_expected_distance") {
        dist <- aldex.expectedDistance(aldex_clr_obj)
        png(filename = opt$output)
        qgraph(dist, layout = "spring", vsize = 1)
        dev.off()
    } else if (opt$analysis_type == "aldex_kw") {
        aldex_obj <- aldex.kw(aldex_clr_obj)
    } else if (opt$analysis_type == "aldex_plot") {
        aldex_obj <- aldex(reads = reads_df,
                           conditions_vector,
                           mc.samples = opt$num_mc_samples,
                           test = opt$aldex_test,
                           effect = effect,
                           include.sample.summary = include_sample_summary,
                           denom = opt$denom,
                           iterate = iterate)
        png(filename = opt$output)
        aldex.plot(x = aldex_obj,
                   type = opt$plot_type,
                   test = opt$plot_test,
                   cutoff.pval = opt$cutoff_pval,
                   cutoff.effect = opt$cutoff_effect,
                   xlab = opt$xlab,
                   ylab = opt$ylab)
        dev.off()
    } else if (opt$analysis_type == "aldex_plot_feature") {
        png(filename = opt$output)
        aldex.plotFeature(aldex_clr_obj, opt$feature_name)
        dev.off()
    } else if (opt$analysis_type == "aldex_ttest") {
        paired_test <- get_boolean_value(opt$paired_test)
        hist_plot <- get_boolean_value(opt$hist_plot)
        aldex_obj <- aldex.ttest(aldex_clr_obj, paired.test = paired_test, hist.plot = hist_plot)
    }
}
if ((opt$analysis_type != "aldex_expected_distance") && (opt$analysis_type != "aldex_plot") && (opt$analysis_type != "aldex_plot_feature")) {
    # Output the ALDEx object.
    write.table(aldex_obj, file = opt$output, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
}
