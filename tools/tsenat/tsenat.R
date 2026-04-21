#!/usr/bin/env Rscript

################################################################################
### TSENAT Galaxy Wrapper Script
### 
### Purpose: Orchestrates the complete TSENAT analysis pipeline
### Inputs: transcript counts, metadata, annotation GFF3
### Outputs: Comprehensive isoform complexity analysis tables and plots
###
### Required Packages: optparse, TSENAT, S4Vectors
################################################################################

suppressPackageStartupMessages({
    library("optparse")
    library("TSENAT")
    library("S4Vectors")
})

# Set locale for consistent output
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Error handling function
error_exit <- function(message, status = 1) {
    cat("ERROR:", message, "\n", file = stderr())
    quit(status = status)
}

# Parse command-line arguments
option_list <- list(
    make_option(c("--salmon_dir"), 
                action = "store", 
                dest = "salmon_dir",
                help = "Path to Salmon output directory (contains sample subdirs with quant.sf files)"),
    make_option(c("--metadata"), 
                action = "store", 
                dest = "metadata",
                help = "Path to sample metadata (TSV with required columns: 'sample', 'condition'; and optional 'subject_col'/'paired_samples' for paired designs)"),
    make_option(c("--annotation"), 
                action = "store", 
                dest = "annotation",
                help = "Path to transcript annotation (GFF3 or GFF3.gz)"),
    make_option(c("--paired"), 
                action = "store_true", 
                dest = "paired",
                default = FALSE,
                help = "Flag indicating paired/repeated-measures design"),
    make_option(c("--control"), 
                action = "store", 
                dest = "control",
                default = "",
                help = "Reference condition level [default: first level]"),
    make_option(c("--q_min"), 
                action = "store", 
                type = "double",
                dest = "q_min",
                default = 0,
                help = "Minimum q value for Tsallis entropy [default: %default]"),
    make_option(c("--q_max"), 
                action = "store", 
                type = "double",
                dest = "q_max",
                default = 2,
                help = "Maximum q value for Tsallis entropy [default: %default]"),
    make_option(c("--q_step"), 
                action = "store", 
                type = "double",
                dest = "q_step",
                default = 0.05,
                help = "Step size for q values [default: %default]"),
    make_option(c("--filter_stringency"), 
                action = "store", 
                dest = "filter_stringency",
                default = "medium",
                help = "Filter stringency: 'soft', 'medium', or 'severe' [default: %default]"),
    
    # Statistical thresholds
    make_option(c("--p_threshold"), 
                action = "store", 
                type = "double",
                dest = "p_threshold",
                default = 0.05,
                help = "P-value significance threshold [default: %default]"),
    make_option(c("--fdr_threshold"), 
                action = "store", 
                type = "double",
                dest = "fdr_threshold",
                default = 0.05,
                help = "FDR significance threshold [default: %default]"),
    make_option(c("--significance_threshold"), 
                action = "store", 
                type = "double",
                dest = "significance_threshold",
                default = 0.05,
                help = "General significance threshold [default: %default]"),
    
    # Bootstrap parameters
    make_option(c("--bootstrap"), 
                action = "store_true", 
                dest = "bootstrap",
                default = FALSE,
                help = "Enable bootstrap confidence intervals"),
    make_option(c("--nboot"), 
                action = "store", 
                type = "integer",
                dest = "nboot",
                default = 1000,
                help = "Number of bootstrap replicates [default: %default]"),
    make_option(c("--bootstrap_method"), 
                action = "store", 
                dest = "bootstrap_method",
                default = "percentile",
                help = "Bootstrap method: 'percentile' or 'bca' [default: %default]"),
    make_option(c("--bootstrap_ci"), 
                action = "store", 
                type = "double",
                dest = "bootstrap_ci",
                default = 0.95,
                help = "Bootstrap confidence level [default: %default]"),
    make_option(c("--bootstrap_include_diagnostics"), 
                action = "store_true", 
                dest = "bootstrap_include_diagnostics",
                default = TRUE,
                help = "Include bootstrap diagnostics [default: %default]"),
    
    # Normalization parameters
    make_option(c("--norm"), 
                action = "store_true", 
                dest = "norm",
                default = TRUE,
                help = "Enable normalization [default: %default]"),
    make_option(c("--norm_method"), 
                action = "store", 
                dest = "norm_method",
                default = NULL,
                help = "Normalization method (NULL for default) [default: %default]"),
    make_option(c("--pseudocount"), 
                action = "store", 
                type = "double",
                dest = "pseudocount",
                default = 0,
                help = "Pseudocount for zeros [default: %default]"),
    make_option(c("--shrinkage"), 
                action = "store", 
                dest = "shrinkage",
                default = "none",
                help = "Shrinkage method: 'none' or other [default: %default]"),
    make_option(c("--min_valid_frac"), 
                action = "store", 
                type = "double",
                dest = "min_valid_frac",
                default = 0.75,
                help = "Minimum valid fraction for filtering [default: %default]"),
    
    # Method selection
    make_option(c("--sait_method"), 
                action = "store", 
                dest = "sait_method",
                default = "gam",
                help = "SAIT fitting method: 'gam', 'lmm', 'fpca', or 'gee' [default: %default]"),
    make_option(c("--sait_pcorr"), 
                action = "store", 
                dest = "sait_pcorr",
                default = "BH",
                help = "P-value correction: 'BH', 'bonferroni', 'hochberg', or 'holm' [default: %default]"),
    make_option(c("--sait_jis_fdr"), 
                action = "store_true", 
                dest = "jis_use_sait_fdr",
                default = TRUE,
                help = "Use SAIT FDR for jackknife isoform switching [default: %default]"),
    make_option(c("--divergence_ci"), 
                action = "store", 
                type = "double",
                dest = "divergence_ci",
                default = 0.95,
                help = "Divergence confidence level [default: %default]"),
    make_option(c("--assumptions_checks"), 
                action = "store", 
                dest = "assumptions_checks",
                default = "all",
                help = "Assumptions checks: 'rank', 'gam', or 'all' [default: %default]"),
    
    # Performance
    make_option(c("--nthreads"), 
                action = "store", 
                type = "integer",
                dest = "nthreads",
                default = 1,
                help = "Number of threads for parallel computation [default: %default]"),
    
    make_option(c("--save_intermediates"), 
                action = "store_true", 
                dest = "save_intermediates",
                default = FALSE,
                help = "Save intermediate analysis objects"),
    make_option(c("--output_prefix"), 
                action = "store", 
                dest = "output_prefix",
                default = "tsenat",
                help = "Prefix for output files [default: %default]"),
    make_option(c("--output_dir"), 
                action = "store", 
                dest = "output_dir",
                default = "./",
                help = "Output directory path [default: %default]"),
    
    # Metadata column names
    make_option(c("--sample_col"), 
                action = "store", 
                dest = "sample_col",
                default = "sample",
                help = "Column name for sample IDs in metadata [default: %default]"),
    make_option(c("--condition_col"), 
                action = "store", 
                dest = "condition_col",
                default = "condition",
                help = "Column name for experimental conditions in metadata [default: %default]"),
    make_option(c("--subject_col"), 
                action = "store", 
                dest = "subject_col",
                default = "",
                help = "Column name for subject/pairing information (for paired designs) [default: auto-detect]"),
    
    # Additional filtering and normalization options
    make_option(c("--min_tpm"), 
                action = "store", 
                type = "double",
                dest = "min_tpm",
                default = NULL,
                help = "Minimum TPM threshold for filtering [default: %default]"),
    make_option(c("--min_isoform_abundance"), 
                action = "store", 
                type = "double",
                dest = "min_isoform_abundance",
                default = NULL,
                help = "Minimum relative isoform abundance within gene [default: %default]"),
    make_option(c("--reference_group"), 
                action = "store", 
                dest = "reference_group",
                default = "",
                help = "Reference group for normalization method 'relative_reference' [default: %default]"),
    make_option(c("--control_group"), 
                action = "store", 
                dest = "control_group",
                default = "",
                help = "Control group for divergence analysis [default: %default]"),
    
    # SAIT advanced options
    make_option(c("--multicorr"), 
                action = "store", 
                dest = "multicorr",
                default = NULL,
                help = "Advanced multiple correction method: 'hochberg', 'westfall-young', 'benjamini-yekutieli' [default: %default]"),
    make_option(c("--corstr"), 
                action = "store", 
                dest = "corstr",
                default = NULL,
                help = "Correlation structure for GEE: 'ar1', 'exchangeable', 'independence' [default: %default]"),
    
    # Subsampling options
    make_option(c("--subset_n_genes"), 
                action = "store", 
                type = "integer",
                dest = "subset_n_genes",
                default = NULL,
                help = "Number of genes to retain (subsampling) [default: %default]"),
    make_option(c("--subset_select_by"), 
                action = "store", 
                dest = "subset_select_by",
                default = "variance",
                help = "Gene selection method:  'variance', 'mean', or 'random' [default: %default]"),
    make_option(c("--subset_seed"), 
                action = "store", 
                type = "integer",
                dest = "subset_seed",
                default = 42,
                help = "Random seed for reproducible gene selection [default: %default]"),
    
    # Diversity output type (advanced)
    make_option(c("--what"), 
                action = "store", 
                dest = "what_output",
                default = "S",
                help = "Diversity output type: 'S' (entropy, default) or 'D' (diversity) [default: %default]"),
    
    # SRH advanced parameters
    make_option(c("--wy_randomizations"), 
                action = "store", 
                type = "integer",
                dest = "wy_randomizations",
                default = 500,
                help = "Westfall-Young randomizations for SRH [default: %default]"),
    make_option(c("--eta2_threshold_moderate"), 
                action = "store", 
                type = "double",
                dest = "eta2_threshold_moderate",
                default = 0.01,
                help = "Moderate effect size threshold for SRH [default: %default]"),
    make_option(c("--eta2_threshold_strong"), 
                action = "store", 
                type = "double",
                dest = "eta2_threshold_strong",
                default = 0.1,
                help = "Strong effect size threshold for SRH [default: %default]"),
    make_option(c("--min_nperm"), 
                action = "store", 
                type = "integer",
                dest = "min_nperm",
                default = 100,
                help = "Minimum permutations for SRH [default: %default]"),
    make_option(c("--max_nperm"), 
                action = "store", 
                type = "integer",
                dest = "max_nperm",
                default = 10000,
                help = "Maximum permutations for SRH [default: %default]"),
    
    # M-estimator parameters
    make_option(c("--loss_type"), 
                action = "store", 
                dest = "loss_type",
                default = "huber",
                help = "M-estimator loss function: 'huber', 'tukey', or 'lsq' [default: %default]"),
    make_option(c("--max_iter"), 
                action = "store", 
                type = "integer",
                dest = "max_iter",
                default = 50,
                help = "M-estimator maximum iterations [default: %default]"),
    make_option(c("--m_tol"), 
                action = "store", 
                type = "double",
                dest = "m_tol",
                default = 1e-6,
                help = "M-estimator convergence tolerance [default: %default]"),
    make_option(c("--q_combine_method"), 
                action = "store", 
                dest = "q_combine_method",
                default = "mean",
                help = "M-estimator q-combination method: 'mean' or 'median' [default: %default]"),
    make_option(c("--influence_threshold"), 
                action = "store", 
                type = "double",
                dest = "influence_threshold",
                default = 0.75,
                help = "M-estimator influence threshold [default: %default]"),
    make_option(c("--scale_method"), 
                action = "store", 
                dest = "scale_method",
                default = "mad",
                help = "M-estimator scale method: 'mad', 'proposal2', or 's-estimator' [default: %default]")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Validate required inputs
if (is.null(args$salmon_dir) || !dir.exists(args$salmon_dir)) {
    error_exit("Required input '--salmon_dir' (Salmon output directory with quant.sf files) not found or not specified")
}
if (is.null(args$metadata) || !file.exists(args$metadata)) {
    error_exit("Required input file 'metadata' not found or not specified")
}
if (is.null(args$annotation) || !file.exists(args$annotation)) {
    error_exit("Required input file 'annotation' not found or not specified")
}

# Validate Salmon directory (build_analysis will handle all loading)
if (!dir.exists(args$salmon_dir)) {
    error_exit(paste("Salmon directory not found:", args$salmon_dir))
}

tryCatch({
    # ============================================================================
    # STEP 1: Load and validate input data
    # ============================================================================
    metadata <- read.table(args$metadata, 
                          header = TRUE, 
                          sep = "\t", 
                          stringsAsFactors = FALSE)
    
    # Validate required columns exist
    if (!args$sample_col %in% colnames(metadata)) {
        stop("Metadata must have '", args$sample_col, "' column (sample column)")
    }
    if (!args$condition_col %in% colnames(metadata)) {
        stop("Metadata must have '", args$condition_col, "' column (condition column)")
    }
    
    # Detect subject_col for paired designs (check for various standard names)
    subject_col_name <- NULL
    if (args$paired) {
        # If user provided a subject_col parameter, use it
        if (args$subject_col != "" && args$subject_col != "NULL") {
            if (args$subject_col %in% colnames(metadata)) {
                subject_col_name <- args$subject_col
            } else {
                stop("Specified subject column '", args$subject_col, "' not found in metadata")
            }
        } else {
            # Auto-detect: check for various standard names
            potential_names <- c("subject_col", "paired_samples", "subject_id", "pair_id", "subject")
            detected <- potential_names[potential_names %in% colnames(metadata)]
            if (length(detected) > 0) {
                subject_col_name <- detected[1]
            } else {
                stop("Paired design specified but no subject/pairing column found (expected: subject_col, paired_samples, subject_id, pair_id, or subject)")
            }
        }
    }
    
    # Validate sample names match
    missing_samples <- setdiff(colnames(readcounts), metadata[[args$sample_col]])
    if (length(missing_samples) > 0) {
        stop("Samples in readcounts not found in metadata: ", 
             paste(missing_samples, collapse = ", "))
    }
    
    # ============================================================================
    # STEP 2: Build TSENAT configuration
    # ============================================================================
    
    # Create q-spectrum
    q_spectrum <- seq(args$q_min, args$q_max, by = args$q_step)
    
    # Set control level if specified
    control <- if (args$control != "" && args$control != "NULL") args$control else NULL
    
    # Build configuration
    config <- TSENAT_config(
        sample_col = args$sample_col,
        condition_col = args$condition_col,
        subject_col = subject_col_name,
        paired = args$paired,
        control = control,
        q = q_spectrum,
        stringency = args$filter_stringency,
        p_threshold = args$p_threshold,
        fdr_threshold = args$fdr_threshold,
        significance_threshold = args$significance_threshold,
        bootstrap = args$bootstrap,
        nboot = args$nboot,
        bootstrap_method = args$bootstrap_method,
        bootstrap_ci = args$bootstrap_ci,
        bootstrap_include_diagnostics = args$bootstrap_include_diagnostics,
        norm = args$norm,
        norm_method = args$norm_method,
        pseudocount = args$pseudocount,
        shrinkage = args$shrinkage,
        min_valid_frac = args$min_valid_frac,
        reference_group = if (args$reference_group != "" && args$reference_group != "NULL") args$reference_group else NULL,
        sait_method = args$sait_method,
        sait_pcorr = args$sait_pcorr,
        multicorr = args$multicorr,
        corstr = args$corstr,
        jis_use_sait_fdr = args$jis_use_sait_fdr,
        divergence_ci = args$divergence_ci,
        control_group = if (args$control_group != "" && args$control_group != "NULL") args$control_group else NULL,
        assumptions_checks = args$assumptions_checks,
        nthreads = args$nthreads,
        subset_n_genes = args$subset_n_genes,
        subset_select_by = args$subset_select_by,
        subset_seed = args$subset_seed,
        min_tpm = args$min_tpm,
        min_isoform_abundance = args$min_isoform_abundance,
        what = args$what_output,
        wy_randomizations = args$wy_randomizations,
        eta2_threshold_moderate = args$eta2_threshold_moderate,
        eta2_threshold_strong = args$eta2_threshold_strong,
        min_nperm = args$min_nperm,
        max_nperm = args$max_nperm,
        loss_type = args$loss_type,
        max_iter = args$max_iter,
        tol = args$m_tol,
        q_combine_method = args$q_combine_method,
        influence_threshold = args$influence_threshold,
        scale_method = args$scale_method
    )
    
    
    # ============================================================================
    # STEP 3: Build analysis object
    # ============================================================================
    
    # Note: salmon_dir is passed directly; build_analysis() will auto-load:
    #   - readcounts (from NumReads column)
    #   - tpm (from TPM column)
    #   - effective_length (from EffectiveLength column)
    analysis <- build_analysis(
        salmon_dir = args$salmon_dir,
        metadata = metadata,
        tx2gene = args$annotation,
        config = config
    )
    
    # ============================================================================
    # STEP 4: Run complete TSENAT pipeline (includes filtering internally)
    # ============================================================================
    
    # Run TSENAT - handles all output generation (tables + plots) automatically
    result <- TSENAT(
        analysis,
        output_dir = args$output_dir,
        save_output = TRUE,
        output_format = "tsv"
    )
    
    # ============================================================================
    # STEP 5: Save analysis object (optional)
    # ============================================================================
    if (args$save_intermediates) {
        saveRDS(result, 
                file = file.path(args$output_dir, 
                                paste0(args$output_prefix, "_analysis_object.rds")))
    }
    
}, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n", file = stderr())
    traceback()
    quit(status = 1)
})

quit(status = 0)
