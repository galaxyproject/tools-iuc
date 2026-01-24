## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----extdata1-----------------------------------------------------------------
sample_id_file <- system.file("extdata/tcga_sample_ids.tsv",
  package = "SplicingFactory"
)

sample_ids <- read.table(sample_id_file)

## ----extdata2-----------------------------------------------------------------
gene_id_file <- system.file("extdata/tcga_gene_ids.tsv",
  package = "SplicingFactory"
)

gene_ids <- read.table(sample_id_file)

## ----setup--------------------------------------------------------------------
library("SplicingFactory")
library("SummarizedExperiment")

# Load dataset
data(tcga_brca_luma_dataset)

# Extract gene names
genes <- tcga_brca_luma_dataset[, 1]

# Extract read count data without gene names
readcounts <- tcga_brca_luma_dataset[, -1]

# Check read count dataset
dim(readcounts)

head(readcounts[, 1:5])

## ----readfilter---------------------------------------------------------------
tokeep <- rowSums(readcounts > 5) > 5

readcounts <- readcounts[tokeep, ]
genes <- genes[tokeep]

## ----laplace------------------------------------------------------------------
laplace_entropy <- calculate_diversity(readcounts, genes,
  method = "laplace",
  norm = TRUE, verbose = TRUE
)

head(assay(laplace_entropy)[, 1:5])

## ----gini---------------------------------------------------------------------
gini_index <- calculate_diversity(readcounts, genes,
  method = "gini",
  verbose = TRUE
)

head(assay(gini_index)[, 1:5])

## ----divplots-----------------------------------------------------------------
library("tidyr")
library("ggplot2")

# Construct data.frame from SummarizedExperiment result
laplace_data <- cbind(assay(laplace_entropy),
  Gene = rowData(laplace_entropy)$genes
)

# Reshape data.frame
laplace_data <- pivot_longer(laplace_data, -Gene,
  names_to = "sample",
  values_to = "entropy"
)

# Add sample type information
laplace_data$sample_type <- apply(
  laplace_data[, 2], 1,
  function(x) {
    ifelse(grepl("_N", x),
      "Normal", "Tumor"
    )
  }
)

# Filter genes with NA entropy values
laplace_data <- drop_na(laplace_data)

# Update gene names and add diversity type column
laplace_data$Gene <- paste0(laplace_data$Gene, "_", laplace_data$sample_type)
laplace_data$diversity <- "Normalized Laplace entropy"

# Construct data.frame from SummarizedExperiment result
gini_data <- cbind(assay(gini_index), Gene = rowData(gini_index)$genes)

# Reshape data.frame
gini_data <- pivot_longer(gini_data, -Gene,
  names_to = "sample",
  values_to = "gini"
)

# Add sample type information
gini_data$sample_type <- apply(
  gini_data[, 2], 1,
  function(x) {
    ifelse(grepl("_N", x),
      "Normal", "Tumor"
    )
  }
)

# Filter genes with NA gini values
gini_data <- drop_na(gini_data)

# Update gene names and add diversity type column
gini_data$Gene <- paste0(gini_data$Gene, "_", gini_data$sample_type)
gini_data$diversity <- "Gini index"

# Plot diversity data
ggplot() +
  geom_density(
    data = laplace_data, alpha = 0.3,
    aes(x = entropy, group = sample, color = diversity)
  ) +
  geom_density(
    data = gini_data, alpha = 0.3,
    aes(x = gini, group = sample, color = diversity)
  ) +
  facet_grid(. ~ sample_type) +
  scale_color_manual(values = c("black", "darkorchid4")) +
  guides(color = FALSE) +
  theme_minimal() +
  labs(
    x = "Diversity values",
    y = "Density"
  )

# Mean entropy calculation across samples for each gene/sample type combination
laplace_entropy_mean <- aggregate(laplace_data$entropy,
  by = list(laplace_data$Gene), mean
)
colnames(laplace_entropy_mean)[2] <- "mean_entropy"
laplace_entropy_mean <- as_tibble(laplace_entropy_mean)

# Add sample type information
laplace_entropy_mean$sample_type <- apply(
  laplace_entropy_mean[, 1], 1,
  function(x) {
    ifelse(grepl("_Normal", x),
      "Normal", "Tumor"
    )
  }
)

# Add diversity type column
laplace_entropy_mean$diversity <- "Normalized Laplace entropy"

# Mean gini calculation across samples for each gene/sample type combination
gini_mean <- aggregate(gini_data$gini, by = list(gini_data$Gene), mean)
colnames(gini_mean)[2] <- "mean_gini"
gini_mean <- as_tibble(gini_mean)

# Add sample type information
gini_mean$sample_type <- apply(
  gini_mean[, 1], 1,
  function(x) {
    ifelse(grepl("_Normal", x),
      "Normal", "Tumor"
    )
  }
)

# Add diversity type column
gini_mean$diversity <- "Gini index"

ggplot() +
  geom_violin(
    data = laplace_entropy_mean, aes(
      x = sample_type, y = mean_entropy,
      fill = diversity
    ),
    alpha = 0.6
  ) +
  geom_violin(
    data = gini_mean, aes(
      x = sample_type, y = mean_gini,
      fill = diversity
    ),
    alpha = 0.6
  ) +
  scale_fill_viridis_d(name = "Diversity") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Samples",
    y = "Diversity"
  )

## -----------------------------------------------------------------------------
# Update the SummarizedExperiment object with a new sample metadata column for
# sample types, as the the object returned by calculate_diversity does not
# contain this information.
colData(laplace_entropy) <- cbind(colData(laplace_entropy),
  sample_type = ifelse(grepl("_N", laplace_entropy$samples),
    "Normal", "Tumor"
  )
)

# Calculate significant entropy changes
entropy_significance <- calculate_difference(
  x = laplace_entropy, samples = "sample_type",
  control = "Normal",
  method = "mean", test = "wilcoxon",
  verbose = TRUE
)

head(entropy_significance)

## -----------------------------------------------------------------------------
entropy_significance$label <- apply(
  entropy_significance[, c(4, 7)], 1,
  function(x) {
    ifelse(abs(x[1]) >= 0.1 & x[2] < 0.05,
      "significant", "non-significant"
    )
  }
)

entropy_significance$mean <- apply(
  entropy_significance[, c(2, 3)], 1,
  function(x) (x[1] + x[2]) / 2
)

ggplot(entropy_significance, aes(x = mean, y = mean_difference)) +
  geom_point(color = "lightgrey", size = 1) +
  geom_point(
    data = entropy_significance[entropy_significance$label == "significant", ],
    color = "red", size = 1
  ) +
  theme_minimal() +
  labs(
    title = "Normalized Laplace entropy",
    subtitle = "Wilcoxon signed rank test",
    x = "Mean entropy",
    y = "Mean difference"
  )

## -----------------------------------------------------------------------------
ggplot(entropy_significance, aes(
  x = mean_difference,
  y = -log10(adjusted_p_values),
  color = label
)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(0.1, -0.1)) +
  theme_minimal() +
  labs(
    title = "Normalized Laplace entropy",
    subtitle = "Wilcoxon signed rank test",
    x = "Mean difference of entropy values",
    y = "-Log10(adjusted p-value"
  )

## -----------------------------------------------------------------------------

sessionInfo()
