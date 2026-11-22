
args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

source(args[1])


### dplot
## This is a comment to test pushing by NPR 07/22/2020
ggplot(dplot$x)
