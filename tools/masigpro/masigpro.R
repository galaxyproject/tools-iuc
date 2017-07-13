#!/usr/bin/env Rscript

# A command-line interface to maSigPro for use with Galaxy
# written by Clemens Blank.
# Thanks to Bjoern Gruening and Michael Love for their DESeq2
# wrapper as a basis to build upon.

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library("maSigPro")
  library("optparse")
  library("mclust")
})

# The following code fixes an error in the stepback function
# of the maSigPro package. This code is hopefully temporary 
# and can be removed if the fix is included in a future 
# version. The stepback function in the maSigPro namespace
# will be overwritten by the following function.
stepback <- function (y = y, d = d, alfa = 0.05, family = gaussian() , epsilon=0.00001) 
{
    lm1 <- glm(y ~ ., data = d, family=family, epsilon=epsilon)
    result <- summary(lm1)
    max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
    if (length(result$coefficients[, 4][-1]) == 1) {
      if (max > alfa) {
        max = 0 
        lm1 <- glm(y ~ 1,  family=family, epsilon=epsilon)
      }
    }
    while (max > alfa) {
        varout <- names(result$coefficients[, 4][-1])[result$coefficients[, 
            4][-1] == max][1]
        pos <- position(matrix = d, vari = varout)
        d <- d[, -pos]
        if (length(result$coefficients[, 4][-1]) == 2) {
            min <- min(result$coefficients[, 4][-1], na.rm = TRUE)
            lastname <- names(result$coefficients[, 4][-1])[result$coefficients[,4][-1] == min]
        }
        if (is.null(dim(d))) {
            d <- as.data.frame(d)
            colnames(d) <- lastname
        }
        lm1 <- glm(y ~ ., data = d, family=family, epsilon=epsilon)
        result <- summary(lm1)
        max <- max(result$coefficients[, 4][-1], na.rm = TRUE)
        if (length(result$coefficients[, 4][-1]) == 1) {
            max <- result$coefficients[, 4][-1]
            if (max > alfa) {
                max = 0
                lm1 <- glm(y ~ 1,  family=family, epsilon=epsilon)
            }
        }
    }
    return(lm1)
}

unlockBinding("stepback", as.environment("package:maSigPro"))
assignInNamespace("stepback", stepback, ns="maSigPro", envir=as.environment("package:maSigPro"))
assign("stepback", stepback, as.environment("package:maSigPro"))
lockBinding("stepback", as.environment("package:maSigPro"))
# End of temporary code to fix stepback.R

options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
 make_option(c("-q", "--quiet"), action="store_false",
 dest="verbose", help="Print little output"),
 make_option(c("-e", "--edesign"), type="character"),
 make_option(c("-d", "--data"), type="character"),
 make_option(c("-o", "--outfile"), type="character"),
 make_option("--degree", type="integer", default=1),
 make_option("--time_col", type="integer", default=1),
 make_option("--repl_col", type="integer", default=2),
 make_option("--qvalue", type="double", default=0.05),
 make_option("--min_obs", type="integer", default=6),
 make_option("--step_method", type="character", default="backward"),
 make_option("--nvar_correction", type="logical", default=FALSE),
 make_option("--alfa", type="double", default=0.05),
 make_option("--rsq", type="double", default=0.7),
 make_option("--vars", type="character", default="groups"),
 make_option("--significant_intercept", type="character", default="dummy"),
 make_option("--cluster_data", type="integer", default=1),
 make_option(c("-k", "--k"), type="integer", default=9),
 make_option("--print_cluster", type="logical", default=FALSE),
 make_option("--cluster_method", type="character", default="hclust"),
 make_option("--distance", type="character", default="cor"),
 make_option("--agglo_method", type="character", default="ward.D"),
 make_option("--iter_max", type="integer", default=500),
 make_option("--color_mode", type="character", default="rainbow"),
 make_option("--show_fit", type="logical", default=TRUE),
 make_option("--show_lines", type="logical", default=TRUE),
 make_option("--cexlab", type="double", default=0.8),
 make_option("--legend", type="logical", default=TRUE)
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

# enforce the following required arguments
if (is.null(opt$edesign)) {
  cat("'edesign' is required\n")
  q(status=1)
}
if (is.null(opt$data)) {
  cat("'data' is required\n")
  q(status=1)
}
if (is.null(opt$outfile)) {
  cat("'outfile' is required\n")
  q(status=1)
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
} else {
  FALSE
}

edesign <- as.matrix(read.table(opt$edesign, header=TRUE, row.names = 1))

data <- read.table(opt$data, header=TRUE, check.names=FALSE)

results <- maSigPro(data, edesign, degree = opt$degree, time.col = opt$time_col,
         repl.col = opt$repl_col, Q = opt$qvalue, min.obs = opt$min_obs,
         step.method = opt$step_method, nvar.correction = opt$nvar_correction,
         alfa = opt$alfa, rsq = opt$rsq, vars = opt$vars,
         significant.intercept = opt$significant_intercept,
         cluster.data = opt$cluster_data, k = opt$k,
         cluster.method = opt$cluster_method, distance = opt$distance,
         agglo.method = opt$agglo_method, iter.max = opt$iter_max,
         color.mode = opt$color_mode, show.fit = opt$show_fit,
         show.lines = opt$show_lines, cexlab = opt$cexlab,
         legend = opt$legend)

if (opt$print_cluster) {
    for (i in 1:length(results$sig.genes)) {
    
    colname <- paste(names(results$sig.genes)[i], "cluster", sep = "_")
    
    results$summary[colname] <- ""
    results$summary[[colname]][1:length(results$sig.genes[[i]]$sig.profiles$`cluster$cut`)] <-
        results$sig.genes[[i]]$sig.profiles$`cluster$cut`
    }
}

filename <- opt$outfile

write.table((results$summary), file=filename, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

cat("Session information:\n\n")

sessionInfo()
