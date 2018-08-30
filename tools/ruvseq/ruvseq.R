# setup R error handling to go to stderr
library("getopt")
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)

setup_cmdline_options <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  spec <- matrix(c(
    "help", "h", 0, "logical",
    "alpha", "a", 1, "double",
    "min_mean_count", "min_c", 1, "double",
    "min_k", "min_k", 1, "double",
    "max_k", "max_k", 1, "double",
    "sample_json", "s", 1, "character",
    "plots" , "p", 1, "character",
    "header", "H", 0, "logical",
    "txtype", "y", 1, "character",
    "tx2gene", "x", 1, "character"), # a space-sep tx-to-gene map or GTF file (auto detect .gtf/.GTF)
    byrow=TRUE, ncol=4)

  opt <- getopt(spec)
  # if help was asked for print a friendly message
  # and exit with a non-zero error code
  if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
  } else {
    load_libraries()
  }
  return(opt)
}

load_libraries <- function() {
  # Allows displaying help without waiting for libraries to load
  library("tools")
  library("jsonlite")
  library("reshape2")
  library("RUVSeq")
  library("RColorBrewer")
  library("tximport")
  library("DESeq2")
  library("ggrepel")
}

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

# Source get_deseq_dataset.R for getting deseq dataset from htseq/featurecounts/tximport
source_local('get_deseq_dataset.R')

# RUVseq function definitions

plot_pca_rle <- function (set, title) {
  x <- pData(set)[,1]
  colors <- brewer.pal(3, "Set2")
  label <- paste0(' for ', title)
  plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
  title(main=paste0("RLE", label))
  plotPCA(set, col=colors[x], cex=1.2)
  title(main=paste0("PCA", label))
}

plot_factors_of_unwanted_variation <- function(set, method, k){
  pd <- pData(set)
  pd['sample'] <- row.names(pd)
  colnames(pd)[1] <- 'condition'
  d = melt(pd, id.vars = c('sample', 'condition'))
  d['x'] <- 1  # There is no information on the X, so we just fake it to be able to do a scatterplot
  print(ggplot(d, aes(x=x, y=value, color=condition, label=sample)) +
  geom_point() +
  ggtitle(paste0('Factors of unwanted variation for method: ', method, ", k=", k)) +
  facet_wrap( ~ variable, scales = "free_x") +
  geom_text_repel() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
  )
}

create_seq_expression_set <- function (dds, min_mean_count) {
  count_values <- counts(dds)
  filter <- apply(count_values, 1, function(x) mean(x) > min_mean_count)
  filtered <- count_values[filter,]
  set = newSeqExpressionSet(as.matrix(count_values),
                            phenoData = data.frame(colData(dds)$condition, row.names=colnames(filtered)))
  plot_pca_rle(set = set, title = 'raw data')
  set <- betweenLaneNormalization(set, which="upper")
  plot_pca_rle(set = set, title = 'upper quartile normalized')
  return(set)
}

get_empirical_control_genes <- function(set, cutoff_p) {
  x <- pData(set)[,1]
  design <- model.matrix(~x, data=pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  top <- topTags(lrt, n=nrow(set))$table
  top_rows <- rownames(top)[which(top$PValue > cutoff_p)]
  empirical <- rownames(set)[which(!(rownames(set) %in% top_rows))]
  return(empirical)
}

ruv_control_gene_method <- function(set, k, control_genes='empirical', cutoff_p=0.2) {
  if (control_genes == 'empirical') {
    control_genes = get_empirical_control_genes(set, cutoff_p=cutoff_p)
  }
  set <- RUVg(set, control_genes, k=k)
  plot_pca_rle(set, paste0("RUVg with empirical control genes, k=", k))
  plot_factors_of_unwanted_variation(set, method="RUVg with empirical control genes", k=k)
  return(set)
}

ruv_residual_method <- function(set, k) {
  genes <- rownames(counts(set))
  x <- pData(set)[,1]
  # Initial edger residuals
  design <- model.matrix(~x, data=pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type="deviance")
  set <- RUVr(set, genes, k=k, res)
  plot_pca_rle(set = set, title = paste0('RUVr using residuals, k=', k))
  plot_factors_of_unwanted_variation(set, method="RUVr using residuals", k=k)
  return(set)
}

ruv_replicate_method <- function (set, k) {
  genes <- rownames(counts(set))
  x <- pData(set)[,1]
  differences <- makeGroups(x)
  set <- RUVs(set, genes, k=k, differences)
  plot_pca_rle(set, paste0('RUVs with replicate samples, k=', k))
  plot_factors_of_unwanted_variation(set, method="RUVs using replicates", k=k)
  return(set)
}

get_differentially_expressed_genes <- function(dds, contrast, alpha=0.01) {
  r <- results(dds, contrast=contrast, alpha=alpha)
  return(rownames(r[which(r$padj < alpha),]))
}

opt <- setup_cmdline_options()
alpha <- opt$alpha
min_k <- opt$min_k
max_k <- opt$max_k
sample_json <- fromJSON(opt$sample_json)
sample_paths <- sample_json$path
sample_names <- sample_json$label
condition <- as.factor(sample_json$condition)
sampleTable <- data.frame(samplename=sample_names,
                          filename = sample_paths,
                          condition=condition)
rownames(sampleTable) <- sample_names

dds <- get_deseq_dataset(sampleTable, header=opt$header, designFormula= ~ condition, tximport=opt$txtype, txtype=opt$txtype, tx2gene=opt$tx2gene)
if (!is.null(opt$plots)) {
  pdf(opt$plots)
}

# Run through the ruvseq variants
set <- create_seq_expression_set(dds, min_mean_count = opt$min_mean_count)
result <- list(no_correction = set)
for (k in seq(min_k, max_k)) {
  result[[paste0('residual_method_k', k)]] <- ruv_residual_method(set, k=k)
  result[[paste0('replicate_method_k', k)]] <- ruv_replicate_method(set, k=k)
  result[[paste0('control_method_k', k)]] <- ruv_control_gene_method(set, k=k, cutoff_p=0.5)
}

for (name in names(result)) {
  if (!startsWith(name, "no_correction")) {
    set <- result[[name]]
    unwanted_variation <- pData(set)
    df <- data.frame(identifier = rownames(unwanted_variation))
    df <- cbind(df, unwanted_variation)
    colnames(df)[2] <- 'condition'
    write.table(df, file=paste0("batch_effects_", name, ".tabular"),  sep="\t", quote=F, row.names=F)
  }
}

# close the plot device
if (!is.null(opt$plots)) {
  cat("closing plot device\n")
  dev.off()
}

cat("Session information:\n\n")
sessionInfo()
