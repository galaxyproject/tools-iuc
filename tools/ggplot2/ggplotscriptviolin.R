# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# We need to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")


# Import library
library("getopt")
library("reshape2")
library("ggplot2")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)


# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
option_specification = matrix(c(
  'input', 'i', 2, 'character',
  'title', 't',2, 'character',
  'ggfill', 's', 2, 'character',
  'ggcolor', 'z', 2, 'character',
  'xlab', 'x', 2, 'character',
  'ylab', 'y', 2, 'character',
  'drawquartiles', 'a', 2, 'character',
  'transform', 'w', 2, 'character',
  'xaxismin', 'e', 2, 'integer',
  'xaxismax', 'f', 2, 'integer',
  'yaxismin', 'g', 2, 'integer',
  'yaxismax', 'h', 2, 'integer',
  'scaling', 'b', 2, 'character',
  'output', 'o', 2, 'character'
  ), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);



# Print options to see what is going on
cat("\n input: ",options$input)
cat("\n title: ",options$title)
#cat("\n size: ",options$size)
cat("\n xlab: ",options$xlab)
cat("\n ylab: ",options$ylab)
cat("\n output: ",options$output)


if(options$scaling == "Automatic"){
    gg_scaley = NULL
} else {
    gg_scaley = ylim(options$yaxismin,options$yaxismax)
    cat("\n yaxismin: ",options$yaxismin)
    cat("\n yaxismax: ",options$yaxismax)
}


integrated <- read.csv(options$input,sep='\t',header=TRUE)
input <- melt(integrated)



if(options$transform == "log2"){
    input["value"] <- log2(input["value"])
}else if(options$transform == "log2plus1"){
    input["value"] <- log2(input["value"]+1)
}else if(options$transform == "log10"){
    input["value"] <- log10(input["value"])
}else if(options$transform == "log10plus1"){
    input["value"] <- log10(input["value"]+1)
    }else{
}

if(options$drawquartiles == "none"){
    gg_quartile = NULL
} else {
    gg_quartile = c(0.25, 0.5, 0.75)
}



ggplot(input,aes(variable,value)) +geom_violin(scale = "area", colour = options$ggcolor, fill = options$ggfill, draw_quantiles = gg_quartile) +
gg_scaley+theme_bw()+xlab(options$xlab)+ylab(options$ylab)+ggtitle(options$title)

ggsave(file=options$output)
