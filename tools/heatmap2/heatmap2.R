# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# We need to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")


# Import library
library("getopt")
library("RColorBrewer")
library("gplots")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)


# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
option_specification = matrix(c(
  'input', 'i', 2, 'character',
  'title', 't', 2, 'character',
  'transform', 'c', 2, 'character',
  'key', 'k', 2, 'character',
  'colorscheme', 'z', 2, 'character',
  'cluster', 'b', 2, 'character',
  'labels', 'a', 2, 'character',
  'scale', 'd', 2, 'character',
  'output', 'o', 2, 'character'
  ), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);



# Print options to see what is going on
cat("\n input: ",options$input)
cat("\n title: ",options$title)
cat("\n output: ",options$output)

input <- read.delim(options$input,sep='\t',header=TRUE)

mat_input <- data.matrix(input[,2:ncol(input)])
rownames(mat_input) <- input[,1]


hclust_fun = function(x) hclust(x, method="complete")
dist_fun = function(x) dist(x, method="maximum")
distfun=dist_fun
hclustfun=hclust_fun
if(options$transform == "none"){
    linput <- mat_input
}else if(options$transform == "log2"){
    linput <- log2(mat_input)
}else if(options$transform == "log2plus1"){
    linput <- log2(mat_input+1)
}else if(options$transform == "log10"){
    linput <- log10(mat_input)
}else if(options$transform == "log10plus1"){
    linput <- log10(mat_input+1)
    }else{
}

if(options$colorscheme == "whrd"){
  colorscale = colfunc <- colorRampPalette(c("white", "red"))
} else if(options$colorscheme == "whblu"){
  colorscale = colfunc <- colorRampPalette(c("white", "blue"))
}else if(options$colorscheme == "blwhre"){
  colorscale = colfunc <- colorRampPalette(c("blue","white", "red"))
}else{
}




if(options$labels== "both"){
  rlabs = NULL
  clabs = NULL
}else if(options$labels== "rows"){
  rlabs = NULL
  clabs = FALSE
}else if(options$labels== "columns"){
  rlabs = FALSE
  clabs = NULL
}else if(options$labels== "none"){
  rlabs = FALSE
  clabs = FALSE
} else{
}




pdf(file="Rplot.pdf")
colorscale

if(options$cluster== "Default"){

heatmap.2(linput,
          distfun=dist_fun, hclustfun=hclust_fun, scale = options$scale, labRow = rlabs, labCol = clabs,
          col=colfunc(50), trace="none", density.info = "none", margins=c(8,2),
          main = options$title, key.xlab= options$key, keysize=1)
}else{
heatmap.2(linput,
          dendrogram="none", scale = options$scale, labRow = rlabs, labCol = clabs,
          col=colfunc(50), trace="none", density.info = "none", margins=c(8,2),
          main = options$title, key.xlab= options$key, keysize=1)
}




dev.off()
