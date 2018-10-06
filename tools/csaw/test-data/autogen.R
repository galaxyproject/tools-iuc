# Generating a random single-end SAM file. This uses the 'simsam.R' script in
# inst/tests, to avoid redundancy. 

source(system.file("tests", "simsam.R", package="csaw"))

set.seed(3894751)
chromos <- c(chrA=1298, chrB=870, chrC=1345)
nreads <- c(1349, 3291)
fnames <- c("rep1", "rep2")

for (curf in 1:length(fnames)) { 
	pos.chr <- sample(length(chromos), nreads[curf], replace=TRUE)
	pos.pos <- str <- integer(nreads[curf])
	for (i in 1:length(chromos)) {
		current<-pos.chr==i
		pos.pos[current]<-round(runif(sum(current), 1, chromos[i]))
		str[current]<-(rbinom(sum(current), 1, 0.5)==1)
	}

	isdup <- rbinom(nreads[curf], 1, 0.8)==0L
	mapq <- round(runif(nreads[curf], 50, 199))
	simsam(fnames[curf], names(chromos)[pos.chr], pos.pos, str, chromos, is.dup=isdup, mapq=mapq)
}

unlink(paste0(fnames, ".sam"))

# I'm not so worried about paired-end data because we've only got one file for that, and I
# don't need to use it so much as an example. Instead, we just sort and compress that file.

petname <- "pet"
samFile <- paste0(petname, ".sam")
tempName <- paste(petname, "_temp", sep="")
tempName <- asBam(samFile, destination=tempName, overwrite=TRUE, indexDestination=FALSE)
newName <- sortBam(tempName, destination=petname)
indexBam(newName)
unlink(tempName)

