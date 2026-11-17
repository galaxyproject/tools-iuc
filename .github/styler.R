#!/usr/bin/env Rscript

library("argparse")
library("styler")

parser <- ArgumentParser(description = "Call styler")
parser$add_argument("dir",
    metavar = "DIR", type = "character",
    help = "File to parse"
)
parser$add_argument("--dry",
    choices = c("off", "on"), default = "on"
)
args <- parser$parse_args()

file_info <- file.info(args$dir)
is_directory <- file_info$isdir

if (is_directory) {
    captured_output <- capture.output({
        result <- style_dir(args$dir, indent_by = 4, dry = args$dry, recursive = TRUE)
    })
} else {
    captured_output <- capture.output({
        result <- style_file(args$dir, indent_by = 4, dry = args$dry)
    })
}

n <- nrow(subset(result, changed == TRUE))
if (n > 0) {
    if (args$dry == "off") {
        print(paste("Changed", n, "files"))
    } else {
        stop(paste("Linting failed for", n, "files"))
    }
}
