library(dada2, quietly = TRUE)
library(ggplot2, quietly = TRUE)

sample_names <- c("F3D0_S188_L001", "F3D141_S207_L001")
fwd <- c("F3D0_S188_L001_R1_001.fastq.gz", "F3D141_S207_L001_R1_001.fastq.gz")
rev <- c("F3D0_S188_L001_R2_001.fastq.gz", "F3D141_S207_L001_R2_001.fastq.gz")

filt_fwd <- c("filterAndTrim_F3D0_R1.fq.gz", "filterAndTrim_F3D141_R1.fq.gz")
filt_rev <- c("filterAndTrim_F3D0_R2.fq.gz", "filterAndTrim_F3D141_R2.fq.gz")

print("filterAndTrim")

for (i in seq_len(fwd)) {
    ftout <- dada2::filterAndTrim(fwd[i], filt_fwd[i], rev[i], filt_rev[i])
    b <- paste(strsplit(fwd[i], ".", fixed = TRUE)[[1]][1], "tab", sep = ".")
    write.table(ftout, b, quote = FALSE, sep = "\t", col.names = NA)
}

# In the test only the 1st data set is used
t <- data.frame()
t <- rbind(t, ftout[1, ])
colnames(t) <- colnames(ftout)
rownames(t) <- rownames(ftout)[1]
write.table(t, "filterAndTrim.tab", quote = FALSE, sep = "\t", col.names = NA)

names(fwd) <- sample_names
names(rev) <- sample_names
names(filt_fwd) <- sample_names
names(filt_rev) <- sample_names

# Plot quality profile (just for one file, Galaxy compares with sim_size)
print("plots")
qp <- dada2::plotQualityProfile(fwd)
ggsave("qualityProfile_fwd.pdf", qp, width = 20, height = 15, units = c("cm"))
qp <- dada2::plotQualityProfile(rev)
ggsave("qualityProfile_rev.pdf", qp, width = 20, height = 15, units = c("cm"))
qp <- dada2::plotQualityProfile(fwd[1])
ggsave("qualityProfile.pdf", qp, width = 20, height = 15, units = c("cm"))

# Plot complexity (just for one file, Galaxy compares with sim_size)

cp <- dada2::plotComplexity(fwd)
ggsave("complexity_fwd.pdf", cp, width = 20, height = 15, units = c("cm"))
cp <- dada2::plotComplexity(rev)
ggsave("complexity_rev.pdf", cp, width = 20, height = 15, units = c("cm"))
cp <- dada2::plotComplexity(fwd[1])
ggsave("complexity.pdf", cp, width = 20, height = 15, units = c("cm"))


# learn Errors
print("learnErrors")
err_fwd <- dada2::learnErrors(filt_fwd)
saveRDS(err_fwd, file = "learnErrors_R1.Rdata")
plot <- dada2::plotErrors(err_fwd)
ggsave("learnErrors_R1.pdf", plot, width = 20, height = 15, units = c("cm"))

err_rev <- dada2::learnErrors(filt_rev)
saveRDS(err_rev, file = "learnErrors_R2.Rdata")
plot <- dada2::plotErrors(err_rev)
ggsave("learnErrors.pdf", plot, width = 20, height = 15, units = c("cm"))

# dada
print("dada")
dada_fwd <- dada2::dada(filt_fwd, err_fwd)
dada_rev <- dada2::dada(filt_rev, err_rev)
for (id in sample_names) {
    saveRDS(dada_fwd[[id]], file = paste("dada_", id, "_R1.Rdata", sep = ""))
    saveRDS(dada_rev[[id]], file = paste("dada_", id, "_R2.Rdata", sep = ""))
}

# merge pairs
print("mergePairs")
merged <- dada2::mergePairs(dada_fwd, filt_fwd, dada_rev, filt_rev)
for (id in sample_names) {
    saveRDS(merged[[id]], file = paste("mergePairs_", id, ".Rdata", sep = ""))
}


# make sequence table
print("makeSequenceTable")
seqtab <- makeSequenceTable(merged)
write.table(t(seqtab), file = "makeSequenceTable.tab", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

reads_per_seqlen <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)
df <- data.frame(length = as.numeric(names(reads_per_seqlen)), count = reads_per_seqlen)
pdf("makeSequenceTable.pdf")
ggplot(data = df, aes(x = length, y = count)) +
    geom_col() +
    theme_bw()
bequiet <- dev.off()

# remove bimera
print("removeBimera")
seqtab_nochim <- dada2::removeBimeraDenovo(seqtab)
write.table(t(seqtab), file = "removeBimeraDenovo.tab", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

# assign taxonomy/species
tl <- "Level1,Level2,Level3,Level4,Level5"
tl <- strsplit(tl, ",")[[1]]

set.seed(42)
print("assignTaxonomyAndSpecies")
taxa <- dada2::assignTaxonomy(seqtab_nochim, "reference.fa.gz", outputBootstraps = TRUE, taxLevels = tl, multithread = 1)

taxa$tax <- dada2::addSpecies(taxa$tax, "reference_species.fa.gz")
write.table(taxa$tax, file = "assignTaxonomyAddspecies.tab", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

write.table(taxa$boot, file = "assignTaxonomyAddspecies_boot.tab", quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)


## Generate extra test data for parameter testing
print("alternatives")
dada2::filterAndTrim(fwd, c("filterAndTrim_single_F3D0_R1.fq.gz", "filterAndTrim_single_F3D141_R1.fq.gz"), rm.phix = TRUE, orient.fwd = "TACGG")

dada2::filterAndTrim(fwd, c("filterAndTrim_single_trimmers_F3D0_R1.fq.gz", "filterAndTrim_single_trimmers_F3D141_R1.fq.gz"), truncQ = 30, truncLen = 2, trimLeft = 150, trimRight = 2)

dada2::filterAndTrim(fwd, c("filterAndTrim_single_filters_F3D0_R1.fq.gz", "filterAndTrim_single_filters_F3D141_R1.fq.gz"), maxLen = 255, minLen = 60, maxN = 100, minQ = 13, maxEE = 1)


merged_nondef <- dada2::mergePairs(dada_fwd, filt_fwd, dada_rev, filt_rev, minOverlap = 8, maxMismatch = 1, justConcatenate = TRUE, trimOverhang = TRUE)
for (id in sample_names) {
    saveRDS(merged_nondef[[id]], file = paste("mergePairs_", id, "_nondefault.Rdata", sep = ""))
}
rb_dada_fwd <- dada2::removeBimeraDenovo(dada_fwd[["F3D0_S188_L001"]])
write.table(rb_dada_fwd, file = "removeBimeraDenovo_F3D0_dada_uniques.tab", quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

rb_merged <- dada2::removeBimeraDenovo(merged, method = "pooled")
saveRDS(rb_merged, file = "removeBimeraDenovo_F3D0_mergepairs.Rdata")

# SeqCounts
get_n <- function(x) {
    sum(dada2::getUniques(x))
}

print("seqCounts ft")
samples <- list()
samples[["F3D0_S188_L001_R1_001.tab"]] <- read.table("F3D0_S188_L001_R1_001.tab", header = TRUE, sep = "\t", row.names = 1)
dname <- "filter"
tdf <- samples[["F3D0_S188_L001_R1_001.tab"]]
names(tdf) <- paste(dname, names(tdf))
tdf <- cbind(data.frame(samples = names(samples)), tdf)
write.table(tdf, "seqCounts_filter.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

samples <- list()
samples[["F3D0_S188_L001_R1_001.tab"]] <- read.table("F3D0_S188_L001_R1_001.tab", header = TRUE, sep = "\t", row.names = 1)
samples[["F3D141_S207_L001_R1_001.tab"]] <- read.table("F3D141_S207_L001_R1_001.tab", header = TRUE, sep = "\t", row.names = 1)
dname <- "filter"
tdf <- samples[["F3D0_S188_L001_R1_001.tab"]]
tdf <- rbind(tdf, samples[["F3D141_S207_L001_R1_001.tab"]])
names(tdf) <- paste(dname, names(tdf))
tdf <- cbind(data.frame(samples = names(samples)), tdf)
write.table(tdf, "seqCounts_filter_both.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("seqCounts dada")
samples <- list()
samples[["dada_F3D0_S188_L001_R1.Rdata"]] <- readRDS("dada_F3D0_S188_L001_R1.Rdata")
samples[["dada_F3D141_S207_L001_R1.Rdata"]] <- readRDS("dada_F3D141_S207_L001_R1.Rdata")
dname <- "dadaF"
tdf <- data.frame(samples = names(samples))
tdf[[dname]] <- sapply(samples, get_n)
write.table(tdf, "seqCounts_dadaF.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("seqCounts mp")
samples <- list()
samples[["mergePairs_F3D0_S188_L001.Rdata"]] <- readRDS("mergePairs_F3D0_S188_L001.Rdata")
samples[["mergePairs_F3D141_S207_L001.Rdata"]] <- readRDS("mergePairs_F3D141_S207_L001.Rdata")
dname <- "merge"
tdf <- data.frame(samples = names(samples))
tdf[[dname]] <- sapply(samples, get_n)
write.table(tdf, "seqCounts_merge.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("seqCounts st")
samples <- list()
samples <- t(as.matrix(read.table("makeSequenceTable.tab", header = TRUE, sep = "\t", row.names = 1)))
dname <- "seqtab"
tdf <- data.frame(samples = row.names(samples))
tdf[[dname]] <- rowSums(samples)
write.table(tdf, "seqCounts_seqtab.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

print("seqCounts rb")
samples <- list()
samples <- t(as.matrix(read.table("removeBimeraDenovo.tab", header = TRUE, sep = "\t", row.names = 1)))
dname <- "nochim"
tdf <- data.frame(samples = row.names(samples))
tdf[[dname]] <- rowSums(samples)
write.table(tdf, "seqCounts_nochim.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
