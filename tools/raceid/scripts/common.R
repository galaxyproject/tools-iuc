# Load libs
arr <- c("methods", "tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit")
res <- lapply(arr, require, character.only = TRUE)


if (is.na(script_dir) || is.na(config_file) || script_dir == "None" || config_file == "None"){
   stop("Could not find script_dir or config_file")
}

source(paste(script_dir, "RaceID_class.R", sep="/"))

# Load params
source(config_file)

# Common functions
message <- function(...){ print(sprintf(...)) }
