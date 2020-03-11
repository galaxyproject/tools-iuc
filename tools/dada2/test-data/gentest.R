library(dada2, quietly=T)
library(ggplot2, quietly=T)

sample.names <- c('F3D0_S188_L001', 'F3D141_S207_L001')
fwd <- c('F3D0_S188_L001_R1_001.fastq.gz', 'F3D141_S207_L001_R1_001.fastq.gz')
rev <- c('F3D0_S188_L001_R2_001.fastq.gz', 'F3D141_S207_L001_R2_001.fastq.gz')

filt.fwd <- c('filterAndTrim_F3D0_R1.fq.gz', 'filterAndTrim_F3D141_R1.fq.gz')
filt.rev <- c('filterAndTrim_F3D0_R2.fq.gz', 'filterAndTrim_F3D141_R2.fq.gz')

print("filterAndTrim")

for(i in 1:length(fwd)){
	ftout <- filterAndTrim(fwd[i], filt.fwd[i], rev[i], filt.rev[i])
    b <- paste(strsplit(fwd[i], ".", fixed=T)[[1]][1], "tab", sep=".")
    write.table(ftout, b, quote=F, sep="\t", col.names=NA)
}

# In the test only the 1st data set is used
t <- data.frame()
t <- rbind(t, ftout[1,])
colnames(t) <- colnames(ftout)
rownames(t) <- rownames(ftout)[1]
write.table(t, "filterAndTrim.tab", quote=F, sep="\t", col.names=NA)

names(fwd) <- sample.names
names(rev) <- sample.names
names(filt.fwd) <- sample.names
names(filt.rev) <- sample.names

# Plot quality profile (just for one file, Galaxy compares with sim_size)
print("plots")
qp <- plotQualityProfile(fwd)
ggsave('qualityProfile_fwd.pdf', qp, width = 20,height = 15,units = c("cm"))
qp <- plotQualityProfile(rev)
ggsave('qualityProfile_rev.pdf', qp, width = 20,height = 15,units = c("cm"))
qp <- plotQualityProfile(fwd[1])
ggsave('qualityProfile.pdf', qp, width = 20,height = 15,units = c("cm"))

# Plot complexity (just for one file, Galaxy compares with sim_size)

cp <- plotComplexity(fwd)
ggsave('complexity_fwd.pdf', cp, width = 20,height = 15,units = c("cm"))
cp <- plotComplexity(rev)
ggsave('complexity_rev.pdf', cp, width = 20,height = 15,units = c("cm"))
cp <- plotComplexity(fwd[1])
ggsave('complexity.pdf', cp, width = 20,height = 15,units = c("cm"))


# learn Errors
print("learnErrors")
err.fwd <- learnErrors(filt.fwd) 
saveRDS(err.fwd, file='learnErrors_R1.Rdata')
plot <- plotErrors(err.fwd)
ggsave('learnErrors_R1.pdf', plot, width = 20,height = 15,units = c("cm"))

err.rev <- learnErrors(filt.rev) 
saveRDS(err.rev, file='learnErrors_R2.Rdata')
plot <- plotErrors(err.rev)
ggsave('learnErrors.pdf', plot, width = 20,height = 15,units = c("cm"))

# dada
print("dada")
dada.fwd <- dada(filt.fwd, err.fwd)
dada.rev <- dada(filt.rev, err.rev)
for( id in sample.names ){
	saveRDS(dada.fwd[[id]], file=paste("dada_", id,"_R1.Rdata", sep=""))
	saveRDS(dada.rev[[id]], file=paste("dada_", id,"_R2.Rdata", sep=""))
}

# merge pairs
print("mergePairs")
merged <- mergePairs(dada.fwd, filt.fwd, dada.rev, filt.rev)
for( id in sample.names ){
	saveRDS(merged[[id]], file=paste("mergePairs_", id,".Rdata", sep=""))
}


# make sequence table
print("makeSequenceTable")
seqtab <- makeSequenceTable(merged)
write.table(t(seqtab), file="makeSequenceTable.tab", quote=F, sep="\t", row.names = T, col.names = NA)

reads.per.seqlen <- tapply(colSums(seqtab), factor(nchar(getSequences(seqtab))), sum)
df <- data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)
pdf( 'makeSequenceTable.pdf' )
ggplot(data=df, aes(x=length, y=count)) +
    geom_col() +
    theme_bw()
bequiet <- dev.off()

# remove bimera
print("removeBimera")
seqtab.nochim <- removeBimeraDenovo(seqtab)
write.table(t(seqtab), file="removeBimeraDenovo.tab", quote=F, sep="\t", row.names = T, col.names = NA)

# assign taxonomy/species
tl <- 'Level1,Level2,Level3,Level4,Level5'
tl <- strsplit(tl, ",")[[1]]

set.seed(42)
print("assignTaxonomyAndSpecies")
taxa <- assignTaxonomy(seqtab.nochim, 'reference.fa.gz', outputBootstraps = T, taxLevels=tl, multithread = 1)

taxa$tax <- addSpecies(taxa$tax, 'reference_species.fa.gz')
write.table(taxa$tax, file = 'assignTaxonomyAddspecies.tab', quote = F, sep = "\t", row.names = T, col.names = NA)

write.table(taxa$boot, file = 'assignTaxonomyAddspecies_boot.tab', quote = F, sep = "\t", row.names = T, col.names = NA)


## Generate extra test data for parameter testing 
print("alternatives")
filterAndTrim(fwd, c('filterAndTrim_single_F3D0_R1.fq.gz', 'filterAndTrim_single_F3D141_R1.fq.gz'), rm.phix = T, orient.fwd = 'TACGG')

filterAndTrim(fwd, c('filterAndTrim_single_trimmers_F3D0_R1.fq.gz', 'filterAndTrim_single_trimmers_F3D141_R1.fq.gz'), truncQ = 30, truncLen = 2, trimLeft = 150, trimRight = 2)

filterAndTrim(fwd, c('filterAndTrim_single_filters_F3D0_R1.fq.gz', 'filterAndTrim_single_filters_F3D141_R1.fq.gz'), maxLen = 255, minLen = 60, maxN = 100, minQ = 13, maxEE = 1)


merged_nondef <- mergePairs(dada.fwd, filt.fwd, dada.rev, filt.rev, minOverlap = 8, maxMismatch = 1, justConcatenate = TRUE, trimOverhang = TRUE)
for( id in sample.names ){
	saveRDS(merged_nondef[[id]], file=paste("mergePairs_", id,"_nondefault.Rdata", sep=""))
}
rb.dada.fwd <- removeBimeraDenovo(dada.fwd[["F3D0_S188_L001"]])
write.table(rb.dada.fwd, file = 'removeBimeraDenovo_F3D0_dada_uniques.tab', quote = F, sep = "\t", row.names = T, col.names = F)

rb.merged <- removeBimeraDenovo(merged, method="pooled")
saveRDS(rb.merged, file='removeBimeraDenovo_F3D0_mergepairs.Rdata')
 
# SeqCounts
getN <- function(x){ sum(getUniques(x)) }

read.uniques <- function ( fname ) {
    p <- read.table(fname, header=F, sep="\t")
    n <-x[,2]
    names(n)<-x[,1]
}


print("seqCounts ft")
samples = list()
samples[["F3D0_S188_L001_R1_001.tab"]] <- read.table("F3D0_S188_L001_R1_001.tab", header=T, sep="\t", row.names=1)
dname <- "filter"
tdf <- samples[["F3D0_S188_L001_R1_001.tab"]]
names(tdf) <- paste( dname, names(tdf) )
tdf <- cbind( data.frame(samples=names( samples )), tdf)
write.table(tdf, "seqCounts_filter.tab", quote=F, sep="\t", row.names = F, col.names = T)

samples = list()
samples[["F3D0_S188_L001_R1_001.tab"]] <- read.table("F3D0_S188_L001_R1_001.tab", header=T, sep="\t", row.names=1)
samples[["F3D141_S207_L001_R1_001.tab"]] <- read.table("F3D141_S207_L001_R1_001.tab", header=T, sep="\t", row.names=1)
dname <- "filter"
tdf <- samples[["F3D0_S188_L001_R1_001.tab"]]
tdf <- rbind(tdf, samples[["F3D141_S207_L001_R1_001.tab"]])
names(tdf) <- paste( dname, names(tdf) )
tdf <- cbind( data.frame(samples=names( samples )), tdf)
write.table(tdf, "seqCounts_filter_both.tab", quote=F, sep="\t", row.names = F, col.names = T)

print("seqCounts dada")
samples = list()
samples[["dada_F3D0_S188_L001_R1.Rdata"]] <- readRDS('dada_F3D0_S188_L001_R1.Rdata')
samples[["dada_F3D141_S207_L001_R1.Rdata"]] <- readRDS('dada_F3D141_S207_L001_R1.Rdata')
dname <- "dadaF"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_dadaF.tab", quote=F, sep="\t", row.names = F, col.names = T)

print("seqCounts mp")
samples = list()
samples[["mergePairs_F3D0_S188_L001.Rdata"]] <- readRDS('mergePairs_F3D0_S188_L001.Rdata')
samples[["mergePairs_F3D141_S207_L001.Rdata"]] <- readRDS('mergePairs_F3D141_S207_L001.Rdata')
dname <- "merge"
tdf <- data.frame( samples = names(samples) )
tdf[[ dname ]] <- sapply(samples, getN)
write.table(tdf, "seqCounts_merge.tab", quote=F, sep="\t", row.names = F, col.names = T)

print("seqCounts st")
samples = list()
samples <- t(as.matrix( read.table("makeSequenceTable.tab", header=T, sep="\t", row.names=1) ))
dname <- "seqtab"
tdf <- data.frame( samples = row.names(samples) )
tdf[[ dname ]] <- rowSums(samples)
write.table(tdf, "seqCounts_seqtab.tab", quote=F, sep="\t", row.names = F, col.names = T)

print("seqCounts rb")
samples = list()
samples <- t(as.matrix( read.table("removeBimeraDenovo.tab", header=T, sep="\t", row.names=1) ))
dname <- "nochim"
tdf <- data.frame( samples = row.names(samples) )
tdf[[ dname ]] <- rowSums(samples)
write.table(tdf, "seqCounts_nochim.tab", quote=F, sep="\t", row.names = F, col.names = T)

