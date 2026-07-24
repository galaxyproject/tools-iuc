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

options(warn = 1)

tryCatch(
    {
        suppressPackageStartupMessages({
            library("optparse")
            library("TSENAT")
            library("S4Vectors")
        })
    },
    error = function(e) {
        cat("ERROR during library loading:\n", file = stderr())
        cat(conditionMessage(e), "\n", file = stderr())
        quit(status = 1)
    }
)

options(readr.show_col_types = FALSE)
options(readr.num_threads = 1)

# Set locale for consistent output
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Error handling function
error_exit <- function(message, status = 1) {
    quit(status = status)
}

# Parse command-line arguments
option_list <- list(
    make_option(c("--salmon_dir"),
        action = "store",
        dest = "salmon_dir",
        default = "",
        help = "Path to directory containing Salmon quantification subdirectories (one per sample, each with quant.sf file). Sample folder names must match the 'sample' column in metadata."
    ),
    make_option(c("--metadata"),
        action = "store",
        dest = "metadata",
        help = "Path to sample metadata (TSV with required columns: 'sample', 'condition'; and optional 'subject_col'/'paired_samples' for paired designs)"
    ),
    make_option(c("--annotation"),
        action = "store",
        dest = "annotation",
        help = "Path to transcript annotation (GFF3 or GFF3.gz)"
    ),
    make_option(c("--paired"),
        action = "store_true",
        dest = "paired",
        default = FALSE,
        help = "Flag indicating paired/repeated-measures design"
    ),
    make_option(c("--control"),
        action = "store",
        dest = "control",
        default = "",
        help = "Reference condition level [default: first level]"
    ),
    make_option(c("--q_min"),
        action = "store",
        type = "double",
        dest = "q_min",
        default = 0,
        help = "Minimum q value for Tsallis entropy [default: %default]"
    ),
    make_option(c("--q_max"),
        action = "store",
        type = "double",
        dest = "q_max",
        default = 2,
        help = "Maximum q value for Tsallis entropy [default: %default]"
    ),
    make_option(c("--q_step"),
        action = "store",
        type = "double",
        dest = "q_step",
        default = 0.05,
        help = "Step size for q values [default: %default]"
    ),
    make_option(c("--filter_stringency"),
        action = "store",
        dest = "filter_stringency",
        default = "medium",
        help = "Filter stringency: 'soft', 'medium', or 'severe' [default: %default]"
    ),

    # Statistical thresholds
    make_option(c("--p_threshold"),
        action = "store",
        type = "double",
        dest = "p_threshold",
        default = 0.05,
        help = "P-value significance threshold [default: %default]"
    ),
    make_option(c("--fdr_threshold"),
        action = "store",
        type = "double",
        dest = "fdr_threshold",
        default = 0.05,
        help = "FDR significance threshold [default: %default]"
    ),
    make_option(c("--significance_threshold"),
        action = "store",
        type = "double",
        dest = "significance_threshold",
        default = 0.05,
        help = "General significance threshold [default: %default]"
    ),

    # Bootstrap parameters
    make_option(c("--bootstrap"),
        action = "store_true",
        dest = "bootstrap",
        default = FALSE,
        help = "Enable bootstrap confidence intervals"
    ),
    make_option(c("--nboot"),
        action = "store",
        type = "integer",
        dest = "nboot",
        default = 1000,
        help = "Number of bootstrap replicates [default: %default]"
    ),
    make_option(c("--bootstrap_method"),
        action = "store",
        dest = "bootstrap_method",
        default = "percentile",
        help = "Bootstrap method: 'percentile' or 'bca' [default: %default]"
    ),
    make_option(c("--bootstrap_ci"),
        action = "store",
        type = "double",
        dest = "bootstrap_ci",
        default = 0.95,
        help = "Bootstrap confidence level [default: %default]"
    ),
    make_option(c("--bootstrap_include_diagnostics"),
        action = "store_true",
        dest = "bootstrap_include_diagnostics",
        default = TRUE,
        help = "Include bootstrap diagnostics [default: %default]"
    ),

    # Normalization parameters
    make_option(c("--norm"),
        action = "store_true",
        dest = "norm",
        default = TRUE,
        help = "Enable normalization [default: %default]"
    ),
    make_option(c("--norm_method"),
        action = "store",
        dest = "norm_method",
        default = NULL,
        help = "Normalization method: NULL (automatic), 'zscore', 'log_odds_ratio', or 'relative_reference'"
    ),
    make_option(c("--pseudocount"),
        action = "store",
        type = "character",
        dest = "pseudocount",
        default = NULL,
        help = "Pseudocount for zeros [default: NULL = use TSENAT default]"
    ),
    make_option(c("--shrinkage"),
        action = "store",
        dest = "shrinkage",
        default = "none",
        help = "Shrinkage method: 'none' or other [default: %default]"
    ),
    make_option(c("--min_valid_frac"),
        action = "store",
        type = "double",
        dest = "min_valid_frac",
        default = 0.75,
        help = "Minimum valid fraction for filtering [default: %default]"
    ),

    # Method selection
    make_option(c("--sait_method"),
        action = "store",
        dest = "sait_method",
        default = "gam",
        help = "SAIT fitting method: 'gam', 'lmm', 'fpca', or 'gee' [default: %default]"
    ),
    make_option(c("--sait_pcorr"),
        action = "store",
        dest = "sait_pcorr",
        default = "BH",
        help = "P-value correction: 'BH', 'bonferroni', 'hochberg', or 'holm' [default: %default]"
    ),
    make_option(c("--sait_jis_fdr"),
        action = "store_true",
        dest = "jis_use_sait_fdr",
        default = TRUE,
        help = "Use SAIT FDR for jackknife isoform switching [default: %default]"
    ),
    make_option(c("--divergence_ci"),
        action = "store",
        type = "double",
        dest = "divergence_ci",
        default = 0.95,
        help = "Divergence confidence level [default: %default]"
    ),
    make_option(c("--assumptions_checks"),
        action = "store",
        dest = "assumptions_checks",
        default = "all",
        help = "Assumptions checks: 'rank', 'gam', or 'all' [default: %default]"
    ),

    # Performance
    make_option(c("--nthreads"),
        action = "store",
        type = "integer",
        dest = "nthreads",
        default = 1,
        help = "Number of threads for parallel computation [default: %default]"
    ),
    make_option(c("--save_intermediates"),
        action = "store_true",
        dest = "save_intermediates",
        default = FALSE,
        help = "Save intermediate analysis objects"
    ),
    make_option(c("--output_prefix"),
        action = "store",
        dest = "output_prefix",
        default = "tsenat",
        help = "Prefix for output files [default: %default]"
    ),
    make_option(c("--output_dir"),
        action = "store",
        dest = "output_dir",
        default = "./",
        help = "Output directory path [default: %default]"
    ),

    # Metadata column names
    make_option(c("--sample_col"),
        action = "store",
        dest = "sample_col",
        default = "sample",
        help = "Column name for sample IDs in metadata [default: %default]"
    ),
    make_option(c("--condition_col"),
        action = "store",
        dest = "condition_col",
        default = "condition",
        help = "Column name for experimental conditions in metadata [default: %default]"
    ),
    make_option(c("--subject_col"),
        action = "store",
        dest = "subject_col",
        default = "",
        help = "Column name for subject/pairing information (for paired designs) [default: auto-detect]"
    ),

    # Additional filtering and normalization options
    make_option(c("--min_tpm"),
        action = "store",
        type = "double",
        dest = "min_tpm",
        default = 1,
        help = "Minimum TPM threshold for filtering [default: %default]"
    ),
    make_option(c("--min_isoform_abundance"),
        action = "store",
        type = "double",
        dest = "min_isoform_abundance",
        default = NULL,
        help = "Minimum relative isoform abundance within gene (NULL for no filtering)"
    ),
    make_option(c("--reference_group"),
        action = "store",
        dest = "reference_group",
        default = "",
        help = "Reference group for normalization method 'relative_reference' [default: %default]"
    ),
    make_option(c("--control_group"),
        action = "store",
        dest = "control_group",
        default = "",
        help = "Control group for divergence analysis [default: %default]"
    ),

    # SAIT advanced options
    make_option(c("--multicorr"),
        action = "store",
        dest = "multicorr",
        default = NULL,
        help = "Advanced multiple correction method: NULL (none), 'hochberg', 'westfall-young', or 'benjamini-yekutieli'"
    ),
    make_option(c("--corstr"),
        action = "store",
        dest = "corstr",
        default = NULL,
        help = "Correlation structure for GEE: NULL (default), 'ar1', 'exchangeable', or 'independence'"
    ),

    # Subsampling options
    make_option(c("--subset_n_genes"),
        action = "store",
        type = "integer",
        dest = "subset_n_genes",
        default = NULL,
        help = "Number of genes to retain for subsampling (NULL for all genes)"
    ),
    make_option(c("--subset_select_by"),
        action = "store",
        dest = "subset_select_by",
        default = "variance",
        help = "Gene selection method:  'variance', 'mean', or 'random' [default: %default]"
    ),
    make_option(c("--subset_seed"),
        action = "store",
        type = "integer",
        dest = "subset_seed",
        default = 42,
        help = "Random seed for reproducible gene selection [default: %default]"
    ),

    # Diversity output type (advanced)
    make_option(c("--what"),
        action = "store",
        dest = "what_output",
        default = "S",
        help = "Diversity output type: 'S' (Tsallis entropy) or 'D' (Hill numbers)"
    ),

    # SRH advanced parameters
    make_option(c("--wy_randomizations"),
        action = "store",
        type = "integer",
        dest = "wy_randomizations",
        default = 500,
        help = "Westfall-Young randomizations for SRH [default: %default]"
    ),
    make_option(c("--eta2_threshold_moderate"),
        action = "store",
        type = "double",
        dest = "eta2_threshold_moderate",
        default = 0.01,
        help = "Moderate effect size threshold for SRH [default: %default]"
    ),
    make_option(c("--eta2_threshold_strong"),
        action = "store",
        type = "double",
        dest = "eta2_threshold_strong",
        default = 0.1,
        help = "Strong effect size threshold for SRH [default: %default]"
    ),
    make_option(c("--min_nperm"),
        action = "store",
        type = "integer",
        dest = "min_nperm",
        default = 100,
        help = "Minimum permutations for SRH [default: %default]"
    ),
    make_option(c("--max_nperm"),
        action = "store",
        type = "integer",
        dest = "max_nperm",
        default = 10000,
        help = "Maximum permutations for SRH [default: %default]"
    ),

    # M-estimator parameters
    make_option(c("--loss_type"),
        action = "store",
        dest = "loss_type",
        default = "huber",
        help = "M-estimator loss function: 'huber', 'tukey', or 'lsq' [default: %default]"
    ),
    make_option(c("--max_iter"),
        action = "store",
        type = "integer",
        dest = "max_iter",
        default = 50,
        help = "M-estimator maximum iterations [default: %default]"
    ),
    make_option(c("--m_tol"),
        action = "store",
        type = "double",
        dest = "m_tol",
        default = 1e-6,
        help = "M-estimator convergence tolerance [default: %default]"
    ),
    make_option(c("--q_combine_method"),
        action = "store",
        dest = "q_combine_method",
        default = "mean",
        help = "M-estimator q-combination method: 'mean' or 'median' [default: %default]"
    ),
    make_option(c("--influence_threshold"),
        action = "store",
        type = "double",
        dest = "influence_threshold",
        default = 0.75,
        help = "M-estimator influence threshold [default: %default]"
    ),
    make_option(c("--scale_method"),
        action = "store",
        dest = "scale_method",
        default = "mad",
        help = "M-estimator scale method: 'mad', 'proposal2', or 's-estimator' [default: %default]"
    ),
    make_option(c("--skip_unmapped"),
        action = "store_true",
        dest = "skip_unmapped",
        default = TRUE,
        help = "Skip unmapped transcripts (not found in annotation). If FALSE and unmapped transcripts exist, analysis will fail. [default: %default]"
    )
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Handle pseudocount: empty/NULL → 0 (default), numbers → as numeric
if (is.null(args$pseudocount) || args$pseudocount == "") {
    args$pseudocount <- 0
} else {
    args$pseudocount <- as.numeric(args$pseudocount)
    if (is.na(args$pseudocount)) {
        error_exit("Pseudocount must be numeric (0-10) or left empty for default")
    }
}

tryCatch(
    {
        # STEP 1: Load and validate input data
        metadata <- read.table(args$metadata,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE
        )

        # Convert column indices (from data_column Galaxy param) to column names
        # data_column returns "1", "2", "3" as strings for column positions
        if (!is.na(suppressWarnings(as.numeric(args$sample_col)))) {
            args$sample_col <- colnames(metadata)[as.numeric(args$sample_col)]
        }
        if (!is.na(suppressWarnings(as.numeric(args$condition_col)))) {
            args$condition_col <- colnames(metadata)[as.numeric(args$condition_col)]
        }
        if (args$subject_col != "" && args$subject_col != "NULL" &&
            !is.na(suppressWarnings(as.numeric(args$subject_col)))) {
            args$subject_col <- colnames(metadata)[as.numeric(args$subject_col)]
        }

        # Set row names to sample column (required by build_analysis())
        rownames(metadata) <- metadata[[args$sample_col]]

        # Detect subject_col for paired designs
        subject_col_name <- NULL
        if (args$paired && args$subject_col != "" && args$subject_col != "NULL") {
            subject_col_name <- args$subject_col
        }

        # STEP 2: Build TSENAT configuration

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

        # Galaxy XML has already properly handled annotation file format detection
        # and created symlink with correct extension (.gz or uncompressed)

        # Set seed for reproducible analysis (matching test workflow)
        set.seed(42)

        analysis <- build_analysis(
            salmon_dir = args$salmon_dir,
            metadata = metadata,
            tx2gene = args$annotation,
            config = config,
            skip = args$skip_unmapped,
            verbose = TRUE
        )

        # STEP 3: Run complete TSENAT pipeline
        # Run TSENAT - handles all output generation (tables + plots) automatically
        result <- TSENAT(
            analysis,
            output_dir = args$output_dir,
            save_output = TRUE,
            output_format = "tsv"
        )

        # STEP 4: Save analysis object (optional)
        if (args$save_intermediates) {
            saveRDS(result,
                file = file.path(
                    args$output_dir,
                    paste0(args$output_prefix, "_analysis_object.rds")
                )
            )
        }
    },
    error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n", file = stderr())
        quit(status = 1)
    }
)

quit(status = 0)
