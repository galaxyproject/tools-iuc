#' Zero-inflated empirical Bayes adjustment of batch effect in compositional
#' feature abundance data
#' 
#' \code{adjust_batch} takes as input a feature-by-sample matrix of microbial
#' abundances, and performs batch effect adjustment given provided batch and
#' optional covariate variables. It returns the batch-adjusted abundance matrix.
#' Additional options and parameters can be passed through the \code{control}
#' parameter as a list (see details).
#' 
#' \code{control} should be provided as a named list of the following components
#' (can be a subset).
#' \describe{
#' \item{zero_inflation}{
#' logical. Indicates whether or not a zero-inflated model should be
#' run. Default to TRUE (zero-inflated model). If set to FALSE then the 
#' correction will be similar to \code{ComBat} as provided in the \code{sva}
#' package.
#' }
#' \item{pseudo_count}{
#' numeric. Pseudo count to add feature_abd before the methods' log
#' transformation. Default to \code{NULL}, in which case \code{adjust_batch} 
#' will set the pseudo count automatically to half of minimal non-zero values in
#' \code{feature_abd}.
#' }
#' \item{diagnostic_plot}{
#' character. Name for the generated diagnostic figure file. Default to 
#' \code{"adjust_batch_diagnostic.pdf"}. Can be set to \code{NULL} in which 
#' case no output will be generated.}
#' \item{conv}{
#' numeric. Convergence threshold for the method's iterative algorithm for 
#' shrinking batch effect parameters. Default to 1e-4.
#' }
#' \item{maxit}{
#' integer. Maximum number of iterations allowed for the method's iterative
#' algorithm. Default to 1000.
#' }
#' \item{verbose}{
#' logical. Indicates whether or not verbose information will be printed.
#' }
#' }
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param covariates name(s) of covariates to adjust for in the batch correction
#' model.
#' @param data data frame of metadata, columns must include batch and covariates
#' (if specified).
#' @param control a named list of additional control parameters. See details.
#'
#' @return a list, with the following components:
#' \describe{
#' \item{feature_abd_adj}{
#' feature-by-sample matrix of batch-adjusted abundances, normalized to the 
#' same per-sample total abundance as feature_abd.
#' }
#' \item{control}{list of additional control parameters used in the function
#' call.
#' }
#' }
#' @export
#' @author Siyuan Ma, \email{siyuanma@@g.harvard.edu}
#' @examples
#' data("CRC_abd", "CRC_meta")
#' CRC_abd_adj <- adjust_batch(feature_abd = CRC_abd, 
#'                             batch = "studyID", 
#'                             covariates = "study_condition",
#'                             data = CRC_meta)$feature_abd_adj
adjust_batch <- function(feature_abd,
                         batch,
                         covariates = NULL,
                         data,
                         control) {
  # Check and construct controls
  control <- match_control(default = control_adjust_batch,
                           control = control)
  verbose <- control$verbose
  
  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
  if(verbose)
    message("feature_abd is ", type_feature_abd)
  # Check metadata data frame
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd,
                           data = data)
  # Check batch and covariates are included in metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  df_batch <- check_metadata(data = data,
                             variables = batch)
  df_covariates <- check_metadata(data = data,
                                  variables = covariates)
  # Check batch variable
  var_batch <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(x = var_batch)
  if(verbose)
    message("Found ", n_batch, " batches")
  
  # Construct batch and covariate model matrices. Check for confounding
  batchmod <- construct_design(data = df_batch, with_intercept = FALSE)
  # For covariates exclude the intercept term because it'ss covered by the
  # batchmod columns
  mod <- construct_design(data = df_covariates, 
                          with_intercept = TRUE)[, -1, drop = FALSE]
  if(!check_rank(design = mod))
    stop("Covariates are confounded!")
  design <- cbind(batchmod, mod)
  if(!check_rank(design = design))
    stop("Covariates and batch are confounded!")
  if(verbose)
    message("Adjusting for ", ncol(mod),
            " covariate(s) or covariate(s) level(s)")
  
  # Transform data for ComBat fit
  if(is.null(control$pseudo_count)) {
    pseudo_count <- set_pseudo(features = feature_abd)
    if(verbose)
      message("Pseudo count is not specified and set to half of minimal ",
              "non-zero value: ",
              format(pseudo_count, digits = 3, scientific = TRUE))
  } else 
    pseudo_count <- check_pseudo_count(control$pseudo_count)
  log_data <- transform_features(
    features = normalize_features(
      features = feature_abd,
      normalization = "TSS",
      pseudo_count = pseudo_count),
    transform = "LOG")
  
  # Identify data to adjust for
  l_ind <- construct_ind(feature_abd = feature_abd,
                         n_batch = n_batch,
                         design = design,
                         zero_inflation = control$zero_inflation)
  batch_no_correct <- apply(!l_ind$ind_gamma, 2, all)
  if(any(batch_no_correct))
    stop(paste0("The following batch(es) either have no present features,",
                " or\nare confounded with the covariates in",
                " features that are present.\nPlease remove them", 
                " from the data before running batch correction:\n",
                paste0(levels(var_batch)[batch_no_correct], collapse = ", ")))
  if(verbose)
    message("Adjusting for (after filtering) ", sum(l_ind$ind_feature),
            " features")
  
  # Standardize data across features
  if(verbose)
    message("Standardizing data across features")
  stand_fit <- fit_stand_feature(s_data = log_data,
                                 design = design,
                                 l_ind = l_ind)
  s_data <- stand_fit$s_data
  l_stand_feature <- stand_fit$l_stand_feature
  
  # Estimate per-batch location and scale parameters
  # and EB hyper-parameters
  if(verbose)
    message("Estimating batch difference parameters and EB priors")
  params_fit <- fit_EB(s_data = s_data, l_stand_feature = l_stand_feature,
                       batchmod = batchmod, n_batch = n_batch,
                       l_ind = l_ind)
  
  # Shrink per-batch location and scale parameters
  if(verbose)
    message("Performing shrinkage adjustments on batch difference parameters")
  params_shrinked <- fit_shrink(s_data = s_data, l_params = params_fit,
                                batchmod = batchmod, n_batch = n_batch,
                                l_ind = l_ind,
                                control = control)
  
  # Adjust the data
  if(verbose)
    message("Performing batch corrections")
  adj_data <- adjust_EB(s_data = s_data, l_params_shrink = params_shrinked,
                        l_stand_feature = l_stand_feature,
                        batchmod = batchmod, n_batch = n_batch,
                        l_ind = l_ind)
  
  # Transform adjusted data back to the original scale
  # For debugging only, this shouldn't happen
  if(any(is.na(adj_data)))
    stop("There are missing values in the adjusted data!")
  feature_abd_adj <- back_transform_abd(adj_data = adj_data,
                                        feature_abd = feature_abd,
                                        type_feature_abd = type_feature_abd)
  
  
  # If required, generate diagnostic plots
  if(!is.null(control$diagnostic_plot))
    diagnostic_adjust_batch(feature_abd = feature_abd,
                            feature_abd_adj = feature_abd_adj,
                            var_batch = var_batch,
                            gamma_hat = params_fit$gamma_hat,
                            gamma_star = params_shrinked$gamma_star,
                            output = control$diagnostic_plot)
  
  return(list(feature_abd_adj = feature_abd_adj,
              control = control))
}