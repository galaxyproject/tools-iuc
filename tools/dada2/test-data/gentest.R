library(dada2, quietly=T)
library(ggplot2, quietly=T)

fwd <- c('F3D0_S188_L001_R1_001.fastq.gz')
rev <- c('F3D0_S188_L001_R2_001.fastq.gz')

sample.names <- c('F3D0_S188_L001')

names(fwd) <- sample.names
names(rev) <- sample.names


filt.fwd <- c('filterAndTrim_F3D0_R1.fq.gz')
filt.rev <- c('filterAndTrim_F3D0_R2.fq.gz')

print("filterAndTrim")

ftout <- filterAndTrim(fwd, filt.fwd, rev, filt.rev)

# In the test no name can be given to the collection
rownames(ftout) <- c( 'Unnamed Collection' )
write.table(ftout, "filterAndTrim_F3D0.tab", quote=F, sep="\t", col.names=NA)

# Plot quality profile (just for one file, Galaxy compares with sim_size)
print("plots")
qp <- plotQualityProfile(fwd)
ggsave('qualityProfile.pdf', qp, width = 20,height = 15,units = c("cm"))

# Plot complexity (just for one file, Galaxy compares with sim_size)

cp <- plotComplexity(fwd)
ggsave('complexity.pdf', cp, width = 20,height = 15,units = c("cm"))


# learn Errors
print("learnErrors")
err.fwd <- learnErrors(filt.fwd) 
saveRDS(err.fwd, file='learnErrors_F3D0_R1.Rdata')
plot <- plotErrors(err.fwd)
ggsave('learnErrors_F3D0_R1.pdf', plot, width = 20,height = 15,units = c("cm"))

err.rev <- learnErrors(filt.fwd) 
saveRDS(err.rev, file='learnErrors_F3D0_R2.Rdata')
plot <- plotErrors(err.rev)
ggsave('learnErrors_F3D0_R2.pdf', plot, width = 20,height = 15,units = c("cm"))

# dada
print("dada")
dada.fwd <- dada(filt.fwd, err.fwd)
saveRDS(dada.fwd, file="dada_F3D0_R1.Rdata")
dada.rev <- dada(filt.rev, err.rev)
saveRDS(dada.rev, file="dada_F3D0_R2.Rdata")

# merge pairs
print("mergePairs")
merged <- mergePairs(dada.fwd, filt.fwd, dada.rev, filt.rev)
saveRDS(merged, file='mergePairs_F3D0.Rdata')

# make sequence table
print("makeSequenceTable")
seqtab <- makeSequenceTable(merged)
write.table(t(seqtab), file="makeSequenceTable_F3D0.tab", quote=F, sep="\t", row.names = T, col.names = NA)

reads.per.seqlen <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)
df <- data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)
pdf( 'makeSequenceTable_F3D0.pdf' )
ggplot(data=df, aes(x=length, y=count)) +
    geom_col() +
    theme_bw()
bequiet <- dev.off()

# remove bimera
print("removeBimera")
seqtab.nochim <- removeBimeraDenovo(seqtab)
write.table(t(seqtab), file="removeBimeraDenovo_F3D0.tab", quote=F, sep="\t", row.names = T, col.names = NA)

# assign taxonomy/species
tl <- 'Level1,Level2,Level3,Level4,Level5'
tl <- strsplit(tl, ",")[[1]]

set.seed(42)
print("assignTaxonomyAndSpecies")
taxa <- assignTaxonomy(seqtab.nochim, 'reference.fa.gz', outputBootstraps = T, taxLevels=c('Level1','Level2','Level3','Level4','Level5'), multithread = 1)

taxa$tax <- addSpecies(taxa$tax, 'reference_species.fa.gz')
write.table(taxa$tax, file = 'assignTaxonomyAddspecies_F3D0.tab', quote = F, sep = "\t", row.names = T, col.names = NA)

write.table(taxa$boot, file = 'assignTaxonomyAddspecies_F3D0_boot.tab', quote = F, sep = "\t", row.names = T, col.names = NA)


## Generate extra test data for parameter testing 
print("alternatives")
filterAndTrim(fwd, c('filterAndTrim_single_F3D0_R1.fq.gz'), rm.phix = T, orient.fwd = 'TACGG')

filterAndTrim(fwd, c('filterAndTrim_single_trimmers_F3D0_R1.fq.gz'), truncQ = 30, truncLen = 2, trimLeft = 150, trimRight = 2)

filterAndTrim(fwd, c('filterAndTrim_single_filters_F3D0_R1.fq.gz'), maxLen = 255, minLen = 60, maxN = 100, minQ = 13, maxEE = 1)


merged_nondef <- mergePairs(dada.fwd, filt.fwd, dada.rev, filt.rev, minOverlap = 8, maxMismatch = 1, justConcatenate = TRUE, trimOverhang = TRUE)
saveRDS(merged_nondef, file='mergePairs_F3D0_nondefault.Rdata')

rb.dada.fwd <- removeBimeraDenovo(dada.fwd)
write.table(rb.dada.fwd, file = 'removeBimeraDenovo_F3D0_dada_uniques.tab', quote = F, sep = "\t", row.names = T, col.names = F)

rb.merged <- removeBimeraDenovo(merged, method="pooled")
saveRDS(rb.merged, file='removeBimeraDenovo_F3D0_mergepairs.Rdata')

# SeqCounts
print("seqCounts")
getN <- function(x){ sum(getUniques(x)) }

read.uniques <- function ( fname ) {
    p <- read.table(fname, header=F, sep="\t")
    n <-x[,2]
    names(n)<-x[,1]
}

samples = list()
samples[["filterAndTrim_F3D0.tab"]] <- read.table("filterAndTrim_F3D0.tab", header=T, sep="\t", row.names=1)
dname <- "filter"
tdf <- samples[["filterAndTrim_F3D0.tab"]]
names(tdf) <- paste( dname, names(tdf) )
tdf <- cbind( data.frame(samples=names( samples )), tdf)
write.table(tdf, "seqCounts_F3D0_filter.tab", quote=F, sep="\t", row.names = F, col.names = T)

samples = list()
samples[["dada_F3D0_R1.Rdata"]] <- readRDS('dada_F3D0_R1.Rdata')
dname <- "dadaF"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_F3D0_dadaF.tab", quote=F, sep="\t", row.names = F, col.names = T)

samples = list()
samples[["mergePairs_F3D0.Rdata"]] <- readRDS('mergePairs_F3D0.Rdata')
dname <- "merge"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_F3D0_merge.tab", quote=F, sep="\t", row.names = F, col.names = T)

samples = list()
samples[["makeSequenceTable_F3D0.tab"]] <- t(as.matrix( read.table("makeSequenceTable_F3D0.tab", header=T, sep="\t", row.names=1) ))
dname <- "seqtab"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_F3D0_seqtab.tab", quote=F, sep="\t", row.names = F, col.names = T)

samples = list()
samples[["removeBimeraDenovo_F3D0.tab"]] <- t(as.matrix( read.table("removeBimeraDenovo_F3D0.tab", header=T, sep="\t", row.names=1) ))
dname <- "nochim"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_F3D0_nochim.tab", quote=F, sep="\t", row.names = F, col.names = T)

