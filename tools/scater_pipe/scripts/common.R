libs <- c("scater", "SingleCellExperiment")
suppressPackageStartupMessages(
    lapply(libs, require, character.only = TRUE)
)

# Intergalactic
if (is.na(script_dir) || is.na(config_file) || script_dir == "None" || config_file == "None"){
   stop("Could not find script_dir or config_file")
}

# Load params
source(config_file)

# Common functions
message <- function(...){ print(sprintf(...)) }
plotter <- function(fname, funct){
    name <- paste(fname, "svg", sep=".")
    svg(name, width=10, height=10)
    force(funct)()
    dev.off()
}

plotterPNG <- function(fname, funct){
    name <- paste(fname, "png", sep=".")
    png(name, width=200, height=200)
    force(funct)()
    dev.off()
}
