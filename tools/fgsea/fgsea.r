options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library("fgsea")
  library("optparse")
})

option_list <- list(
  make_option(c("-rnk_file", "--rnk_file"), type="character", help="Path to file with differential gene expression result"),
  make_option(c("-out_tab","--out_tab"), type="character", help="Path to output file."),
  make_option(c("-gmt_file", "--gmt_file"), type="character", help = "Path to Broad gmt file"),
  make_option(c("-min_size", "--min_size"), type="integer", help="Minimal size of a gene set to test. All pathways below the threshold are excluded."),
  make_option(c("-max_size", "--max_size"), type="integer", help="Maximal size of a gene set to test. All pathways above the threshold are excluded."),
  make_option(c("n_perm", "--n_perm"), type="integer", help="Number of permutations to do. Minimial possible nominal p-value is about 1/nperm"),
  make_option(c("rdaOpt", "--rdaOpt"), type="logical", help="Output RData file"),
  make_option(c("plotOpt", "--plotOpt"), type="logical", help="Output Table plot")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
rnk_file = args$rnk_file
gmt_file = args$gmt_file
out_tab = args$out_tab
min_size = args$min_size
max_size = args$max_size
n_perm = args$n_perm
rdaOpt = args$rdaOpt
plotOpt = args$plotOpt
## Basically using the steps from the fgsea vignette

rankTab <- read.table(rnk_file, header=TRUE, colClasses = c("character", "numeric"))

ranks <-rankTab[,2]
names(ranks) <- rankTab[,1]

## Report an error if gmt_file not found
if(file.exists(gmt_file)) {
  pathways <- gmtPathways(gmt_file)
} else {
  cat(paste("Could not find file ", gmt_file))
  cat("Printing contents of directory:")
  list.files()
}


fgseaRes <- fgsea(pathways, ranks, minSize=min_size, maxSize=max_size, nperm=n_perm)

# Convert leadingEdge column from list to character to output
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, toString)

write.table(fgseaRes, out_tab, sep="\t", row.names=FALSE)

# Make table plot for top pathways

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

## Output Gsea Table plot
if (plotOpt) {
  pdf("GseaTable.pdf")
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
  dev.off()
}

## Output RData file
if (rdaOpt) {
    save.image(file = "fgsea_analysis.RData")
}


