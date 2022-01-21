options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library("goseq")
  library("optparse")
  library("dplyr")
  library("ggplot2")
})

sessionInfo()

option_list <- list(
  make_option("--dge_file", type = "character", help = "Path to file with differential gene expression result"),
  make_option("--length_file", type = "character", default = NULL, help = "Path to tabular file mapping gene id to length"),
  make_option("--genome", type = "character", default = NULL, help = "Genome [used for looking up correct gene length]"),
  make_option("--gene_id", type = "character", default = NULL, help = "Gene ID format of genes in DGE file"),
  make_option("--fetch_cats", type = "character", default = NULL, help = "Categories to get can include one or more of GO:CC, GO:BP, GO:MF, KEGG"),
  make_option("--category_file", type = "character", default = NULL, help = "Path to tabular file with gene_id <-> category mapping"),
  make_option("--wallenius_tab", type = "character", default = NULL, help = "Path to output file with P-values estimated using wallenius distribution"),
  make_option("--nobias_tab", type = "character", default = NULL, help = "Path to output file with P-values estimated using hypergeometric distribution and no correction for gene length bias"),
  make_option("--repcnt", type = "integer", default = 0, help = "Number of repeats for sampling"),
  make_option("--sampling_tab", type = "character", default = NULL, help = "Path to output file with P-values estimated using sampling distribution"),
  make_option("--p_adj_method", type = "character", default = "BH", help = "Multiple hypothesis testing correction method to use"),
  make_option("--use_genes_without_cat", type = "logical", default = FALSE, help = "A large number of gene may have no GO term annotated. If this option is set to FALSE, genes without category will be ignored in the calculation of p-values(default behaviour). If TRUE these genes will count towards the total number of genes outside the tested category (default behaviour prior to version 1.15.2)."),
  make_option("--top_plot", type = "character", default = NULL, help = "Path to output PDF with top10 over-rep GO terms"),
  make_option("--make_plots", default = FALSE, type = "logical", help = "Produce diagnostic plots?"),
  make_option("--length_bias_plot", type = "character", default = NULL, help = "Path to length-bias plot"),
  make_option("--sample_vs_wallenius_plot", type = "character", default = NULL, help = "Path to plot comparing sampling with wallenius p-values"),
  make_option("--rdata", type = "character", default = NULL, help = "Path to RData output file"),
  make_option("--categories_genes_out_fp", type = "character", default = NULL, help = "Path to file with categories (GO/KEGG terms) and associated DE genes")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

if (!is.null(args$fetch_cats)) {
  fetch_cats <- unlist(strsplit(args$fetch_cats, ","))
} else {
  fetch_cats <- "Custom"
}

# format DE genes into named vector suitable for goseq
# check if header is present
first_line <- read.delim(args$dge_file, header = FALSE, nrow = 1)
second_col <- toupper(first_line[, ncol(first_line)])
if (second_col == TRUE || second_col == FALSE) {
  dge_table <- read.delim(args$dge_file, header = FALSE, sep = "\t")
} else {
  dge_table <- read.delim(args$dge_file, header = TRUE, sep = "\t")
}
genes <- as.numeric(as.logical(dge_table[, ncol(dge_table)])) # Last column contains TRUE/FALSE
names(genes) <- dge_table[, 1] # Assuming first column contains gene names

# gene lengths, assuming last column
first_line <- read.delim(args$length_file, header = FALSE, nrow = 1)
if (is.numeric(first_line[, ncol(first_line)])) {
  length_table <- read.delim(args$length_file, header = FALSE, sep = "\t", check.names = FALSE)
} else {
  length_table <- read.delim(args$length_file, header = TRUE, sep = "\t", check.names = FALSE)
}
row.names(length_table) <- length_table[, 1]
# get vector of gene length in same order as the genes
gene_lengths <- length_table[names(genes), ][, ncol(length_table)]

# Estimate PWF
if (args$make_plots) {
  pdf(args$length_bias_plot)
}
pwf <- nullp(genes, genome = args$genome, id = args$gene_id, bias.data = gene_lengths, plot.fit = args$make_plots)
if (args$make_plots) {
  dev.off()
}

# Fetch GO annotations if category_file hasn't been supplied:
if (is.null(args$category_file)) {
  go_map <- getgo(genes = names(genes), genome = args$genome, id = args$gene_id, fetch.cats = fetch_cats)
} else {
  # check for header: first entry in first column must be present in genes, else it's a header
  first_line <- read.delim(args$category_file, header = FALSE, nrow = 1)
  if (first_line[, 1] %in% names(genes)) {
    go_map <- read.delim(args$category_file, header = FALSE)
  } else {
    go_map <- read.delim(args$category_file, header = TRUE)
  }
}

results <- list()

run_goseq <- function(pwf, genome, gene_id, goseq_method, use_genes_without_cat, repcnt, gene2cat, p_adj_method, out_fp) {
  out <- goseq(pwf, genome = genome, id = gene_id, method = goseq_method, use_genes_without_cat = use_genes_without_cat, gene2cat = go_map)
  out$p_adjust_over_represented <- p.adjust(out$over_represented_pvalue, method = p_adj_method)
  out$p_adjust_under_represented <- p.adjust(out$under_represented_pvalue, method = p_adj_method)
  write.table(out, out_fp, sep = "\t", row.names = FALSE, quote = FALSE)
  return(out)
}

# wallenius approximation of p-values
if (!is.null(args$wallenius_tab)) {
  results[["Wallenius"]] <- run_goseq(
    pwf,
    genome = args$genome,
    gene_id = args$gene_id,
    goseq_method = "Wallenius",
    use_genes_without_cat = args$use_genes_without_cat,
    repcnt = args$repcnt,
    gene2cat = go_map,
    p_adj_method = args$p_adj_method,
    out_fp = args$wallenius_tab
  )
}


# hypergeometric (no length bias correction)
if (!is.null(args$nobias_tab)) {
  results[["Hypergeometric"]] <- run_goseq(
    pwf,
    genome = args$genome,
    gene_id = args$gene_id,
    goseq_method = "Hypergeometric",
    use_genes_without_cat = args$use_genes_without_cat,
    repcnt = args$repcnt,
    gene2cat = go_map,
    p_adj_method = args$p_adj_method,
    out_fp = args$nobias_tab
  )
}

# Sampling distribution
if (args$repcnt > 0) {
  results[["Sampling"]] <- run_goseq(
    pwf,
    genome = args$genome,
    gene_id = args$gene_id,
    goseq_method = "Sampling",
    use_genes_without_cat = args$use_genes_without_cat,
    repcnt = args$repcnt,
    gene2cat = go_map,
    p_adj_method = args$p_adj_method,
    out_fp = args$sampling_tab
  )

  # Compare sampling with wallenius
  if (args$make_plots & !is.null(args$wallenius_tab)) {
    pdf(args$sample_vs_wallenius_plot)
    plot(log10(results[["Wallenius"]][, 2]),
      log10(results[["Sampling"]][match(results[["Sampling"]][, 1], results[["Wallenius"]][, 1]), 2]),
      xlab = "log10(Wallenius p-values)",
      ylab = "log10(Sampling p-values)",
      xlim = c(-3, 0)
    )
    abline(0, 1, col = 3, lty = 2)
    dev.off()
  }
}

# Plot the top 10
if (!is.null(args$top_plot)) {
  cats_title <- gsub("GO:", "", args$fetch_cats)
  # modified from https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  pdf(args$top_plot)
  for (m in names(results)) {
    p <- results[[m]] %>%
      top_n(10, wt = -over_represented_pvalue) %>%
      mutate(hitsPerc = numDEInCat * 100 / numInCat) %>%
      ggplot(aes(
        x = hitsPerc,
        y = reorder(substr(term, 1, 40), -over_represented_pvalue), # only use 1st 40 chars of terms otherwise squashes plot
        colour = p_adjust_over_represented,
        size = numDEInCat
      )) +
      geom_point() +
      expand_limits(x = 0) +
      labs(x = "% DE in category", y = "Category", colour = "Adj P value", size = "Count", title = paste("Top over-represented categories in", cats_title), subtitle = paste(m, " method")) +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    print(p)
  }
  dev.off()
}

# Extract the genes to the categories (GO/KEGG terms)
if (!is.null(args$categories_genes_out_fp)) {
  cat2gene <- split(rep(names(go_map), sapply(go_map, length)), unlist(go_map, use.names = FALSE))
  # extract categories (GO/KEGG terms) for all results
  categories <- c()
  for (m in names(results)) {
    categories <- c(categories, results[[m]]$category)
  }
  categories <- unique(categories)
  # extract the DE genes for each catge term
  categories_genes <- data.frame(category = categories, de_genes = rep("", length(categories)))
  categories_genes$de_genes <- as.character(categories_genes$de_genes)
  rownames(categories_genes) <- categories
  for (cat in categories) {
    tmp <- pwf[cat2gene[[cat]], ]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    categories_genes[cat, "de_genes"] <- paste(tmp, collapse = ",")
  }
  # output
  write.table(categories_genes, args$categories_genes_out_fp, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Output RData file
if (!is.null(args$rdata)) {
  save.image(file = args$rdata)
}
