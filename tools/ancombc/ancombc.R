#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ANCOMBC"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--phyloseq"), action = "store", dest = "phyloseq", help = "File containing a phyloseq object"),
    make_option(c("--formula"), action = "store", dest = "formula", help = "Formula"),
    make_option(c("--p_adj_method"), action = "store", dest = "p_adj_method", help = "Method to adjust p-values"),
    make_option(c("--zero_cut"), action = "store", dest = "zero_cut", type = "double", help = "Minimum taxa prevalence"),
    make_option(c("--lib_cut"), action = "store", dest = "lib_cut", type = "integer", help = "Thrshold for filtering samples based on library sizes"),
    make_option(c("--group"), action = "store", dest = "group", help = "Name of the group variable in the metadata"),
    make_option(c("--struc_zero"), action = "store", dest = "struc_zero", help = "Detect structural zeros based on group"),
    make_option(c("--neg_lb"), action = "store", dest = "neg_lb", help = "Classify a taxon as a structural zero using its asymptotic lower bound"),
    make_option(c("--tol"), action = "store", dest = "tol", type = "double", help = "Iteration convergence tolerance for the E-M algorithm"),
    make_option(c("--max_iter"), action = "store", dest = "max_iter", help = "Maximum number of iterations for the E-M algorithm"),
    make_option(c("--conserve"), action = "store", dest = "conserve", help = "Use a conservative variance estimator for the test statistic"),
    make_option(c("--alpha"), action = "store", dest = "alpha", help = "Level of significance"),
    make_option(c("--global"), action = "store", dest = "global", help = "Perform global test"),
    make_option(c("--output_dir"), action = "store", dest = "output_dir", help = "Output directory")
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

get_file_path <- function(dir, file_name) {
    file_path <- paste(dir, file_name, sep = "/")
    return(file_path)
}

write_data_frame <- function(dir, file_name, data_frame) {
    file_path <- get_file_path(dir, file_name)
    write.table(data_frame, file = file_path, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
}

# Convert boolean values to boolean.
struc_zero <- get_boolean_value(opt$struc_zero)
neg_lb <- get_boolean_value(opt$neg_lb)
conserve <- get_boolean_value(opt$conserve)
global <- get_boolean_value(opt$global)

# Construct a phyloseq object.
phyloseq_obj <- readRDS(opt$phyloseq)

# Construct an ANCOM-BC object.
ancombc_obj <- ancombc(
    phyloseq = phyloseq_obj,
    formula = opt$formula,
    p_adj_method = opt$p_adj_method,
    zero_cut = opt$zero_cut,
    lib_cut = opt$lib_cut,
    group = opt$group,
    struc_zero = struc_zero,
    neg_lb = neg_lb,
    tol = opt$tol,
    max_iter = opt$max_iter,
    conserve = conserve,
    alpha = opt$alpha,
    global = global
)

res <- ancombc_obj$res

# Write the outputs.
write_data_frame(opt$output_dir, "feature_table.tabular", ancombc_obj$feature_table)
write_data_frame(opt$output_dir, "zero_ind.tabular", ancombc_obj$zero_ind)
write.csv2(ancombc_obj$samp_frac, file = get_file_path(opt$output_dir, "samp_frac.tabular"), row.names = FALSE, col.names = FALSE, sep = "\t")
write_data_frame(opt$output_dir, "resid.tabular", ancombc_obj$resid)
write(ancombc_obj$delta_em, file = get_file_path(opt$output_dir, "delta_em.tabular"))
write(ancombc_obj$delta_wls, file = get_file_path(opt$output_dir, "delta_wls.tabular"))
write_data_frame(opt$output_dir, "res_beta.tabular", res$beta)
write_data_frame(opt$output_dir, "res_se.tabular", res$se)
write_data_frame(opt$output_dir, "res_W.tabular", res$W)
write_data_frame(opt$output_dir, "res_p_val.tabular", res$p_val)
write_data_frame(opt$output_dir, "res_q_val.tabular", res$q_val)
write_data_frame(opt$output_dir, "res_diff_abn.tabular", res$diff_abn)
write_data_frame(opt$output_dir, "res_global.tabular", ancombc_obj$res_global)
