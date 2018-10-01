#
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(dplyr)
  library(getopt)
  library(ggplot2)
  library(ggrepel)
})

options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

spec <- matrix(c(
  "input", "i", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Below from http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html

results = read.table(opt$input, header=TRUE, sep="\t")
results = mutate(results, sig=ifelse(results$padj<0.05, "FDR<0.05", "Not Sig"))

pdf("out.pdf")
p = ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black")) +
  geom_text_repel(data=filter(results, padj<0.05), aes(label=Gene))
print(p)
dev.off()

cat("Session information:\n\n")
sessionInfo()

