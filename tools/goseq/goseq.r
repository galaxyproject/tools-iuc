options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("goseq")
    library("optparse")
    library("dplyr")
    library("ggplot2")
})

option_list <- list(
    make_option(c("-d", "--dge_file"), type="character", help="Path to file with differential gene expression result"),
    make_option(c("-w","--wallenius_tab"), type="character", help="Path to output file with P-values estimated using wallenius distribution."),
    make_option(c("-s","--sampling_tab"), type="character", default=FALSE, help="Path to output file with P-values estimated using sampling distribution."),
    make_option(c("-n","--nobias_tab"), type="character", default=FALSE, help="Path to output file with P-values estimated using hypergeometric distribution and no correction for gene length bias."),
    make_option(c("-l","--length_bias_plot"), type="character", default=FALSE, help="Path to length-bias plot."),
    make_option(c("-sw","--sample_vs_wallenius_plot"), type="character", default=FALSE, help="Path to plot comparing sampling with wallenius p-values."),
    make_option(c("-r", "--repcnt"), type="integer", default=100, help="Number of repeats for sampling"),
    make_option(c("-lf", "--length_file"), type="character", default="FALSE", help = "Path to tabular file mapping gene id to length"),
    make_option(c("-cat_file", "--category_file"), default="FALSE", type="character", help = "Path to tabular file with gene_id <-> category mapping."),
    make_option(c("-g", "--genome"), default=NULL, type="character", help = "Genome [used for looking up correct gene length]"),
    make_option(c("-i", "--gene_id"), default=NULL, type="character", help = "Gene ID format of genes in DGE file"),
    make_option(c("-p", "--p_adj_method"), default="BH", type="character", help="Multiple hypothesis testing correction method to use"),
    make_option(c("-cat", "--use_genes_without_cat"), default=FALSE, type="logical",
                help="A large number of gene may have no GO term annotated. If this option is set to FALSE, genes without category will be ignored in the calculation of p-values(default behaviour). If TRUE these genes will count towards the total number of genes outside the tested category (default behaviour prior to version 1.15.2)."),
    make_option(c("-plots", "--make_plots"), default=FALSE, type="logical", help="produce diagnostic plots?"),
    make_option(c("-fc", "--fetch_cats"), default=NULL, type="character", help="Categories to get can include one or more of GO:CC, GO:BP, GO:MF, KEGG"),
    make_option(c("-rd", "--rdata"), default=NULL, type="character", help="Path to RData output file."),
    make_option(c("-tp", "--top_plot"), default=NULL, type="logical", help="Output PDF with top10 over-rep GO terms?")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
dge_file = args$dge_file
category_file = args$category_file
length_file = args$length_file
genome = args$genome
gene_id = args$gene_id
wallenius_tab = args$wallenius_tab
sampling_tab = args$sampling_tab
nobias_tab = args$nobias_tab
length_bias_plot = args$length_bias_plot
sample_vs_wallenius_plot = args$sample_vs_wallenius_plot
repcnt = args$repcnt
p_adj_method = args$p_adj_method
use_genes_without_cat = args$use_genes_without_cat
make_plots = args$make_plots
rdata = args$rdata

if (!is.null(args$fetch_cats)) {
  fetch_cats = unlist(strsplit(args$fetch_cats, ","))
} else {
  fetch_cats = "Custom"
}

# format DE genes into named vector suitable for goseq
# check if header is present
first_line = read.delim(dge_file, header = FALSE, nrow=1)
second_col = toupper(first_line[, ncol(first_line)])
if (second_col == TRUE || second_col == FALSE) {
    dge_table = read.delim(dge_file, header = FALSE, sep="\t")
} else {
    dge_table = read.delim(dge_file, header = TRUE, sep="\t")
}
genes = as.numeric(as.logical(dge_table[,ncol(dge_table)])) # Last column contains TRUE/FALSE
names(genes) = dge_table[,1] # Assuming first column contains gene names

# gene lengths, assuming last column
if (length_file != "FALSE" ) {
  first_line = read.delim(length_file, header = FALSE, nrow=1)
  if (is.numeric(first_line[, ncol(first_line)])) {
    length_table = read.delim(length_file, header=FALSE, sep="\t", check.names=FALSE)
    } else {
    length_table = read.delim(length_file, header=TRUE, sep="\t", check.names=FALSE)
    }
  row.names(length_table) = length_table[,1]
  gene_lengths = length_table[names(genes),][,ncol(length_table)]
  } else {
  gene_lengths = getlength(names(genes), genome, gene_id)
  }

# Estimate PWF

if (make_plots != 'false') {
  pdf(length_bias_plot)
}
pwf=nullp(genes, genome = genome, id = gene_id, bias.data = gene_lengths, plot.fit=make_plots)
if (make_plots != 'false') {
  dev.off()
}

# Fetch GO annotations if category_file hasn't been supplied:
if (category_file == "FALSE") {
  go_map=getgo(genes = names(genes), genome=genome, id=gene_id, fetch.cats=fetch_cats)
  } else {
  # check for header: first entry in first column must be present in genes, else it's a header
  first_line = read.delim(category_file, header = FALSE, nrow=1)
  if (first_line[,1] %in% names(genes)) {
     go_map = read.delim(category_file, header = FALSE)
     } else {
     go_map = read.delim(category_file, header= TRUE)
    }
}

results <- list()

# wallenius approximation of p-values
if (wallenius_tab != FALSE) {
  GO.wall=goseq(pwf, genome = genome, id = gene_id, use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)
  GO.wall$p.adjust.over_represented = p.adjust(GO.wall$over_represented_pvalue, method=p_adj_method)
  GO.wall$p.adjust.under_represented = p.adjust(GO.wall$under_represented_pvalue, method=p_adj_method)
  write.table(GO.wall, args$wallenius_tab, sep="\t", row.names = FALSE, quote = FALSE)
  results[['Wallenius']] <- GO.wall
}

# hypergeometric (no length bias correction)
if (nobias_tab != FALSE) {
  GO.nobias=goseq(pwf, genome = genome, id = gene_id, method="Hypergeometric", use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)
  GO.nobias$p.adjust.over_represented = p.adjust(GO.nobias$over_represented_pvalue, method=p_adj_method)
  GO.nobias$p.adjust.under_represented = p.adjust(GO.nobias$under_represented_pvalue, method=p_adj_method)
  write.table(GO.nobias, args$nobias_tab, sep="\t", row.names = FALSE, quote = FALSE)
  results[['Hypergeometric']] <- GO.nobias
}

# Sampling distribution
if (repcnt > 0) {

  # capture the sampling progress so it doesn't fill stdout  
  zz <- file("/dev/null", open = "wt")
  sink(zz)
  GO.samp=goseq(pwf, genome = genome, id = gene_id, method="Sampling", repcnt=repcnt, use_genes_without_cat = use_genes_without_cat, gene2cat=go_map)
  sink()
  
  GO.samp$p.adjust.over_represented = p.adjust(GO.samp$over_represented_pvalue, method=p_adj_method)
  GO.samp$p.adjust.under_represented = p.adjust(GO.samp$under_represented_pvalue, method=p_adj_method)
  write.table(GO.samp, sampling_tab, sep="\t", row.names = FALSE, quote = FALSE)
  # Compare sampling with wallenius
  if (make_plots == TRUE) {
  pdf(sample_vs_wallenius_plot)
  plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
     abline(0,1,col=3,lty=2)
  dev.off()
  }
  results[['Sampling']] <- GO.samp
}

if (!is.null(args$top_plot)) {
  cats_title <- gsub("GO:","", args$fetch_cats)
  # modified from https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  pdf("top10.pdf")
  for (m in names(results)) {
    p <- results[[m]] %>%
      top_n(10, wt=-p.adjust.over_represented)  %>%
      mutate(hitsPerc=numDEInCat*100/numInCat) %>%
      ggplot(aes(x=hitsPerc,
                   y=substr(term, 1, 40), # only use 1st 40 chars of terms otherwise squashes plot
                   colour=p.adjust.over_represented,
                   size=numDEInCat)) +
      geom_point() +
      expand_limits(x=0) +
      labs(x="% DE in category", y="Category", colour="adj. P value", size="Count", title=paste("Top over-represented categories in", cats_title), subtitle=paste(m, " method")) +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    print(p)
  }
  dev.off()
}

# Output RData file
if (!is.null(args$rdata)) {
  save.image(file = "goseq_analysis.RData")
}


sessionInfo()
