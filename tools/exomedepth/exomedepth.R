# Load ExomeDepth library (without warnings)
suppressMessages(library(ExomeDepth))

# Import parameters from xml wrapper (args_file)
args  <- commandArgs(trailingOnly=TRUE)
param <- read.table(args[1],sep="=", as.is=TRUE)

# Set common parameters
target      <- param[match("target",param[,1]),2]
trans_prob  <- as.numeric(param[match("trans_prob",param[,1]),2])
output      <- param[match("output",param[,1]),2]
test_vs_ref <- as.logical(param[match("test_vs_ref",param[,1]),2])

# Create symbolic links for multiple bam and bai 
bam         <- param[param[,1]=="bam",2]
bam_bai     <- param[param[,1]=="bam_bai",2]
bam_label   <- param[param[,1]=="bam_label",2]
bam_label   <- gsub(" ", "_", bam_label)

for(i in 1:length(bam)){
  stopifnot(file.symlink(bam[i], paste(bam_label[i], "bam", sep=".")))
  stopifnot(file.symlink(bam_bai[i], paste(bam_label[i], "bam.bai", sep=".")))
}

# Generate read count data
BAMFiles <- paste(bam_label, "bam", sep=".")
sink("/dev/null")
ExomeCount <- suppressMessages(getBamCounts(bed.file=target, bam.files = BAMFiles))
sink()

# Convert counts in a data frame
ExomeCount.dafr <- as(ExomeCount[, colnames(ExomeCount)], 'data.frame')

# Prepare the main matrix of read count data
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern='.bam')])

# Remove .bam from sample name
colnames(ExomeCount.mat) <- gsub(".bam", "", colnames(ExomeCount.mat))

# Set nsamples == 1 if mode is test vs reference, assuming test is sample 1
nsamples <- ifelse(test_vs_ref, 1, ncol(ExomeCount.mat))

# Loop over samples
for (i in 1:nsamples){

	# Create the aggregate reference set for this sample
	my.choice <- suppressWarnings(suppressMessages(
					select.reference.set(test.counts = ExomeCount.mat[,i],  
                           				 reference.counts = subset(ExomeCount.mat, select=-i),  
                           			     bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                           		         n.bins.reduced = 10000)))

	my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop=FALSE],
                               	   MAR = 1,
                               	   FUN = sum)
                               
	# Now create the ExomeDepth object for the CNVs call
	all.exons<-suppressWarnings(suppressMessages(new('ExomeDepth',
               test = ExomeCount.mat[,i],
               reference = my.reference.selected,
               formula = 'cbind(test,reference)~1')))


	# Now call the CNVs
	result <- try(all.exons<-suppressMessages(CallCNVs(x=all.exons,
                            transition.probability = trans_prob,
                            chromosome = ExomeCount.dafr$space,
                            start = ExomeCount.dafr$start,
                            end = ExomeCount.dafr$end,
                            name = ExomeCount.dafr$names)), silent=T)

	# Next if CNVs are not detected
	if (class(result)=="try-error"){
		next
	}

	# Compute correlation between ref and test
	my.cor <- cor(all.exons@reference, all.exons@test)
	n.call <- nrow(all.exons@CNV.calls)

	# Write results
	my.results <- cbind(all.exons@CNV.calls[,c(7,5,6,3)], 
            		 	sample=colnames(ExomeCount.mat)[i],
            			corr=my.cor,
            			all.exons@CNV.calls[,c(4,9,12)])
            			
    	# Re-order by chr and position
    	chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
    	my.results[,1] <- factor(my.results[,1], levels=chrOrder)
    	my.results <- my.results[order(my.results[,1], my.results[,2], my.results[,3]),]
    
	write.table(sep='\t', quote=FALSE, file = output,
            		x = my.results,
       			row.names = FALSE, col.names = FALSE, dec=".", append=TRUE)
}
