library(tidyverse)
library(ggupset)

args <- commandArgs(trailingOnly = TRUE)
n_int <- as.integer(args[1])
x_lab <- as.character(args[2])
y_lab <- as.character(args[3])
width <- as.integer(args[4])
height <- as.integer(args[5])
files <- tail(args, -5)

# gene, (list of matching conditions)
data <- list()
for (idx in seq_along(files)) {
    k <- files[idx]
    data[k] <- read_tsv(k, col_names = c("genes"))
}

all_genes <- c()
for (gene_list in data) {
    for (gene in gene_list) {
        all_genes <- append(all_genes, gene)
    }
}
all_genes <- unique(sort(all_genes))

trans <- tibble(gene = character(), files = character())
for (gene_idx in seq_along(all_genes)) {
    gene <- all_genes[gene_idx]
    files <- names(data)[grep(gene, data)]
    files <- paste(files, sep = ",", collapse = ",")
    trans <- trans %>% add_row(gene = gene, files = files)
}

write_tsv(trans, "upset.tsv")
trans$f <- str_split(trans$files, ",")

pl <- trans %>% ggplot(aes(x = f)) +
    geom_bar() +
    scale_x_upset(n_intersections = n_int, order_by = "freq") +
    xlab(x_lab) +
    ylab(y_lab)
ggsave("upset-freq.png", width = width, height = height, units = "px")
pl <- trans %>% ggplot(aes(x = f)) +
    geom_bar() +
    scale_x_upset(n_intersections = n_int, order_by = "degree") +
    xlab(x_lab) +
    ylab(y_lab)
ggsave("upset-degree.png", width = width, height = height, units = "px")
