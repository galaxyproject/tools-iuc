# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# We need to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")


# Import library
library(getopt)
library(Rtsne)
library(ggplot2)
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)


# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
option_specification = matrix(c(
  'input', 'i', 2, 'character',
  'title', 't', 2, 'character',
  'transform', 'c', 2, 'character',
  'perp', 'p', 2, 'double',
  'seed', 's', 2, 'double',
  'data', 'd', 2, 'double',
  'name', 'n', 2, 'double',
  'theta', 'z', 2, 'double',
  'legend', 'l', 2, 'character',
  'output', 'o', 2, 'character'
  ), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);



# Print options to see what is going on
cat("\n input: ",options$input)
cat("\n title: ",options$title)
cat("\n name: ",options$name)
cat("\n data: ",options$data)
cat("\n seed: ",options$seed)
cat("\n perp: ",options$perp)
cat("\n theta: ",options$theta)
cat("\n output: ",options$output)

all <- read.delim(options$input,sep='\t',header=TRUE)

data = options$data
name = options$name

if(options$transform == "none"){
    log.all <- all[, data:ncol(all)]
}else if(options$transform == "log2"){
log.all <- log2(all[, data:ncol(all)])
}else if(options$transform == "log2plus1"){
log.all <- log2(1+all[, data:ncol(all)])
}else if(options$transform == "log10"){
log.all <- log(all[, data:ncol(all)])
}else if(options$transform == "log10plus1"){
log.all <- log(1+all[, data:ncol(all)])
    }else{
}

#Show/hide legend
if(options$legend == "yes"){
    gg_legend = NULL
} else {
    gg_legend = theme(legend.position="none")
}

set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(as.matrix(log.all[,1:ncol(log.all)]), perplexity = options$perp, theta = options$theta) # Run TSNE

embedding <- as.data.frame(tsne_out$Y)
embedding$Class <- as.factor(sub("Class_", "", all[,name]))

ggplot(embedding, aes(x=V1, y=V2, color=Class)) +
     geom_point(size=1.25) + gg_legend +
     xlab("") + ylab("") +
     ggtitle(options$title)

ggsave("Rplot.pdf")