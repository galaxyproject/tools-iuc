# Load ExomeDepth library (without warnings)
suppressMessages(library(ExomeDepth))

# Import parameters from xml wrapper (args_file)
args <- commandArgs(trailingOnly = TRUE)
param <- read.table(args[1], sep = "=", as.is = TRUE)

# Set common parameters
target <- param[match("target", param[, 1]), 2]
trans_prob <- as.numeric(param[match("trans_prob", param[, 1]), 2])
output <- param[match("output", param[, 1]), 2]
test_vs_ref <- as.logical(param[match("test_vs_ref", param[, 1]), 2])

# Create symbolic links for multiple bam and bai
bam <- param[param[, 1] == "bam", 2]
bam_bai <- param[param[, 1] == "bam_bai", 2]
bam_label <- param[param[, 1] == "bam_label", 2]
bam_label <- gsub(" ", "_", bam_label)

for (i in seq_len(length(bam))) {
    stopifnot(file.symlink(bam[i], paste(bam_label[i], "bam", sep = ".")))
    stopifnot(file.symlink(bam_bai[i], paste(bam_label[i], "bam.bai", sep = ".")))
}

# Generate read count data
bam_files <- paste(bam_label, "bam", sep = ".")
sink("/dev/null")
exome_count <- suppressMessages(getBamCounts(bed.file = target, bam.files = bam_files))
sink()

# Convert counts in a data frame
exome_count_dafr <- as(exome_count[, colnames(exome_count)], "data.frame")

# Prepare the main matrix of read count data
exome_count_mat <- as.matrix(exome_count_dafr[, grep(names(exome_count_dafr), pattern = ".bam")])

# Remove .bam from sample name
colnames(exome_count_mat) <- gsub(".bam", "", colnames(exome_count_mat))

# Set nsamples == 1 if mode is test vs reference, assuming test is sample 1
nsamples <- ifelse(test_vs_ref, 1, ncol(exome_count_mat))

# Loop over samples
for (i in 1:nsamples) {
    # Create the aggregate reference set for this sample
    my_choice <- suppressWarnings(suppressMessages(select.reference.set(
        test.counts = exome_count_mat[, i],
        reference.counts = subset(exome_count_mat, select = -i),
        bin.length = (exome_count_dafr$end - exome_count_dafr$start) / 1000,
        n.bins.reduced = 10000
    )))

    my_reference_selected <- apply(
        X = exome_count_mat[, my_choice$reference.choice, drop = FALSE],
        MAR = 1,
        FUN = sum
    )

    # Now create the ExomeDepth object for the CNVs call
    all.exons <- suppressWarnings(suppressMessages(
        new("ExomeDepth",
        test = exome_count_mat[, i],
        reference = my_reference_selected,
        formula = "cbind(test,reference)~1"
    )))

    # Now call the CNVs
    result <- try(all.exons <- suppressMessages(CallCNVs(
        x = all.exons,
        transition.probability = trans_prob,
        chromosome = exome_count_dafr$space,
        start = exome_count_dafr$start,
        end = exome_count_dafr$end,
        name = exome_count_dafr$names
    )), silent = T)

    # Next if CNVs are not detected
    if (class(result) == "try-error") {
        next
    }

    # Compute correlation between ref and test
    my_cor <- cor(all.exons@reference, all.exons@test)

    # Write results
    my_results <- cbind(
        all.exons@CNV.calls[, c(7, 5, 6, 3)],
        sample = colnames(exome_count_mat)[i],
        corr = my_cor,
        all.exons@CNV.calls[, c(4, 9, 12)]
    )

    # Re-order by chr and position
    chr_order <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
    my_results[, 1] <- factor(my_results[, 1], levels = chr_order)
    my_results <- my_results[order(my_results[, 1], my_results[, 2], my_results[, 3]), ]

    write.table(
        sep = "\t",
        quote = FALSE,
        file = output,
        x = my_results,
        row.names = FALSE,
        col.names = FALSE,
        dec = ".",
        append = TRUE
    )
}
