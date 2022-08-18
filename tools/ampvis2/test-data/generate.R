library(ampvis2)

# subset taxa using 200 random OTUs
aalborgwwtps <- amp_subset_taxa(AalborgWWTPs, tax_vector = sample(AalborgWWTPs$tax$OTU, 200))

ape::write.tree(aalborgwwtps$tree, "AalborgWWTPs.nwk")
amp_export_fasta(aalborgwwtps, "AalborgWWTPs.fa", tax = FALSE)
amp_export_otutable(aalborgwwtps, "AalborgWWTPs.otu")
write.table(aalborgwwtps$tax, file = "AalborgWWTPs.tax", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(aalborgwwtps$metadata, file = "AalborgWWTPs.tsv", quote = FALSE, sep = "\t")

# construct data for merge_ampvis2
t <- amp_load(otutable = "AalborgWWTPs.otu.csv", metadata = "AalborgWWTPs.tsv", taxonomy = "AalborgWWTPs.tax.tsv", fasta = "AalborgWWTPs.fa", tree = "AalborgWWTPs.nwk")
d_2010 <- amp_subset_samples(t, Period %in% "2010")
# 60 samples and 109 OTUs have been filtered
# Before: 67 samples and 200 OTUs
# After: 7 samples and 91 OTUs
d_2011 <- amp_subset_samples(t, Period %in% "2011")
# 61 samples and 94 OTUs have been filtered
# Before: 67 samples and 200 OTUs
# After: 6 samples and 106 OTUs

saveRDS(d_2010, "AalborgWWTPs.2010.rds")
saveRDS(d_2011, "AalborgWWTPs.2011.rds")

# => merging should give 13 samples
