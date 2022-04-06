library(ampvis2)

# subset taxa using 200 random OTUs
aalborgwwtps <- amp_subset_taxa(AalborgWWTPs, tax_vector = sample(AalborgWWTPs$tax$OTU, 200))

ape::write.tree(aalborgwwtps$tree, "AalborgWWTPs.nwk")
amp_export_fasta(aalborgwwtps, "AalborgWWTPs.fa", tax = F)
amp_export_otutable(aalborgwwtps, "AalborgWWTPs.otu")
write.table(aalborgwwtps$tax, file = "AalborgWWTPs.tax", quote = F, sep = "\t", row.names = F)
write.table(aalborgwwtps$metadata, file = "AalborgWWTPs.tsv", quote = F, sep = "\t")
