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
  'size', 's', 2, 'double',
  'xlab', 'x', 2, 'character',
  'ylab', 'y', 2, 'character',
  'binwidth', 'b', 2, 'double',
  'xaxismin', 'e', 2, 'integer',
  'xaxismax', 'f', 2, 'integer',
  'yaxismin', 'g', 2, 'integer',
  'yaxismax', 'h', 2, 'integer',
  'colors', 'q', 2, 'character',
  'colorscheme', 'z', 2, 'character',
  'colororder', 'r', 2, 'integer',
  'facet', 'a', 2, 'character',
  'transform', 'c', 2, 'character',
  'woutputdim', 'w', 2, 'integer',
  'houtputdim', 'd', 2, 'integer',
  'dim', 'k', 2, 'character',
  'scaling', 'j', 2, 'character',
  'legend', 'l', 2, 'character',
  'density', 'm', 2, 'character',
  'output', 'o', 2, 'character'
  ), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);



# Print options to see what is going on
cat("\n input: ",options$input)
cat("\n title: ",options$title)
cat("\n size: ",options$size)
cat("\n transform: ",options$alpha)
cat("\n xlab: ",options$xlab)
cat("\n xlab: ",options$xlab)
cat("\n binwidth: ",options$binwidth)
cat("\n output: ",options$output)

integrated <- read.csv(options$input,sep='\t',header=TRUE)
input <- melt(integrated)

#Show/hide legend
if(options$legend == "yes"){
    gg_legend = NULL
} else {
    gg_legend = theme(legend.position="none")
}

#density 
if(options$density == "counts"){
    gg_density = ggplot(input,aes(value, color=variable))
    gg_freq = NULL
}else if(options$density == "nfreq"){
    gg_density = ggplot(input,aes(value, ..ncount.., color=variable))
    gg_freq = NULL
}else if(options$density == "freq"){
    gg_density = ggplot(input,aes(value, color=variable))
    gg_freq = aes(y=..count../sum(..count..))
}else {

}

#Choose between automatically scaled x and y axis or user defined
if(options$scaling == "Automatic"){
    gg_scalex = NULL
    gg_scaley = NULL
} else {
    gg_scalex = xlim(options$xaxismin,options$xaxismax)
    gg_scaley = ylim(options$yaxismin,options$yaxismax)
    cat("\n xaxismin: ",options$xaxismin)
    cat("\n xaxismax: ",options$xaxismax)
    cat("\n yaxismin: ",options$yaxismin)
    cat("\n yaxismax: ",options$yaxismax)
}

#Choose dimensions of output pdf
if(options$dim == "Default"){
    gg_width = 7
    gg_height = 7
} else {
    gg_width = options$woutputdim
    gg_height = options$houtputdim 
}


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


if(options$facet == "facet"){
    gg_facet = facet_wrap(~ variable)
  } else {
    gg_facet = NULL
  }

if(options$colorscheme == "Default"){
    gg_colorscale = NULL
  } else {
    gg_colorscale = scale_color_brewer(palette=options$colors, direction=options$colororder)
  }

gg_density +
geom_freqpoly(gg_freq,binwidth=options$binwidth,size=options$size)+gg_facet+gg_colorscale+
gg_scalex+gg_scaley+theme_bw()+xlab(options$xlab)+ylab(options$ylab)+gg_legend+
ggtitle(options$title)

ggsave(file=options$output)
