options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-rnk_file", "--rnk_file"), type="character", help="Path to ranked genes file"),
  make_option(c("-header", "--header"), type="logical", help = "Does ranked genes file have a header"),
  make_option(c("-sets_file", "--sets_file"), type="character", help = "Path to gene sets file"),
  make_option(c("-gmt", "--gmt"), type="logical", help = "Is the sets file in GMT format"),
  make_option(c("-out_tab","--out_tab"), type="character", help="Path to output file"),
  make_option(c("-min_size", "--min_size"), type="integer", help="Minimal size of a gene set to test. All pathways below the threshold are excluded."),
  make_option(c("-max_size", "--max_size"), type="integer", help="Maximal size of a gene set to test. All pathways above the threshold are excluded."),
  make_option(c("-n_perm", "--n_perm"), type="integer", help="Number of permutations to do. Minimial possible nominal p-value is about 1/nperm"),
  make_option(c("-rda_opt", "--rda_opt"), type="logical", help="Output RData file"),
  make_option(c("-plot_opt", "--plot_opt"), type="logical", help="Output plot"),
  make_option(c("-top_num", "--top_num"), type="integer", help="Top number of pathways to plot")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
rnk_file = args$rnk_file
if (args$header) {
  header = TRUE
} else {
  header = FALSE
}
sets_file = args$sets_file
gmt = args$gmt
out_tab = args$out_tab
min_size = args$min_size
max_size = args$max_size
n_perm = args$n_perm
rda_opt = args$rda_opt
plot_opt = args$plot_opt
top_num = args$top_num

## Basically using the steps from the fgsea vignette
rankTab <- read.table(rnk_file, header=header, colClasses = c("character", "numeric"))

ranks <-rankTab[,2]
names(ranks) <- rankTab[,1]

if (gmt) {
  pathways <- gmtPathways(sets_file)
} else {
  pathways <- load(sets_file)
  pathways <- get(pathways)
}

# need to set seed for reproducibility https://github.com/ctlab/fgsea/issues/12
set.seed(42)
fgseaRes <- fgsea(pathways, ranks, minSize=min_size, maxSize=max_size, nperm=n_perm)
fgseaRes <- fgseaRes[order(pval), ]
# Convert leadingEdge column from list to character to output
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, toString)

write.table(fgseaRes, out_tab, sep="\t", row.names=FALSE, quote=FALSE)

if (plot_opt) {
  pdf("fgsea_plots.pdf", width=8)

  topPathways <- head(fgseaRes, n=top_num)
  topPathways <- topPathways$pathway

  ## Make summary table plot for top pathways
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5,
  colwidths = c(5.3,3,0.7, 0.9, 0.9))

  # Make enrichment plots for top pathways
  for (i in topPathways) {
    p <- plotEnrichment(pathways[[i]], ranks) + labs(title=i)
    print(p)
  }

  dev.off()
}

## Output RData file
if (rda_opt) {
    save.image(file = "fgsea_analysis.RData")
}