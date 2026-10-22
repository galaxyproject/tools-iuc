suppressPackageStartupMessages(library(microDecon))
suppressPackageStartupMessages(library(optparse))

# Define command-line options
option_list <- list(
    make_option(c("-m", "--mode"), type = "character", help = "Mode of operation: decon, remove.cont, remove.thresh, decon.diff", metavar = "MODE"),
    make_option(c("-d", "--data_file"), type = "character", help = "Path to the data file (CSV format expected)", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output file from remove.cont or remove.thresh (used in decon.diff)", metavar = "FILE"),
    make_option(c("-b", "--numb_blanks"), type = "integer", default = NULL, help = "Number of blank samples"),
    make_option(c("-n", "--numb_ind"), type = "character", default = NULL, help = "Number of individuals (eval-parsed)"),
    make_option(c("-t", "--taxa"), type = "logical", default = FALSE, help = "Taxa flag (TRUE or FALSE)"),
    make_option(c("-r", "--runs"), type = "integer", default = NULL, help = "Number of runs"),
    make_option(c("-T", "--thresh"), type = "double", default = NULL, help = "Threshold value"),
    make_option(c("-p", "--prop_thresh"), type = "double", default = NULL, help = "Proportional threshold value"),
    make_option(c("-g", "--regression"), type = "double", default = NULL, help = "Regression value"),
    make_option(c("-l", "--low_threshold"), type = "double", default = 40, help = "Low threshold [default: %default]"),
    make_option(c("-u", "--up_threshold"), type = "double", default = 400, help = "Upper threshold [default: %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Read main input data
microbe_data <- read.csv(opt$data_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (opt$mode == "decon") {
    result <- decon(
        data = microbe_data,
        numb.blanks = opt$numb_blanks,
        numb.ind = eval(parse(text = opt$numb_ind)),
        taxa = opt$taxa,
        runs = opt$runs,
        thresh = opt$thresh,
        prop.thresh = opt$prop_thresh,
        regression = opt$regression,
        low.threshold = opt$low_threshold,
        up.threshold = opt$up_threshold
    )
    write.csv(result$decon.table, "decon_table.csv", row.names = FALSE)
    write.csv(result$reads.removed, "reads_removed.csv", row.names = FALSE)
    write.csv(result$sum.per.group, "difference_sum.csv", row.names = FALSE)
    write.csv(result$mean.per.group, "difference_mean.csv", row.names = FALSE)
    write.csv(result$OTUs.removed, "OTUs_removed.csv", row.names = FALSE)
} else if (opt$mode == "remove_cont") {
    result <- remove.cont(
        data = microbe_data,
        numb.blanks = opt$numb_blanks,
        taxa = opt$taxa,
        runs = opt$runs,
        regression = opt$regression,
        low.threshold = opt$low_threshold,
        up.threshold = opt$up_threshold
    )
    write.csv(result, "decon_table.csv", row.names = FALSE)
} else if (opt$mode == "remove_thresh") {
    result <- remove.thresh(
        data = microbe_data,
        numb.ind = eval(parse(text = opt$numb_ind)),
        taxa = opt$taxa,
        thresh = opt$thresh,
        prop.thresh = opt$prop_thresh
    )
    write.csv(result, "decon_table.csv", row.names = FALSE)
} else if (opt$mode == "decon_diff") {
    if (is.null(opt$output)) stop("Error: --output must be provided for decon.diff mode")
    output_data <- read.csv(opt$output, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    result <- decon.diff(
        data = microbe_data,
        output = output_data,
        numb.blanks = opt$numb_blanks,
        numb.ind = eval(parse(text = opt$numb_ind)),
        taxa = opt$taxa
    )
    write.csv(result$decon.table, "decon_table.csv", row.names = FALSE)
    write.csv(result$reads.removed, "reads_removed.csv", row.names = FALSE)
    write.csv(result$sum.per.group, "difference_sum.csv", row.names = FALSE)
    write.csv(result$mean.per.group, "difference_mean.csv", row.names = FALSE)
    write.csv(result$OTUs.removed, "OTUs_removed.csv", row.names = FALSE)
}
