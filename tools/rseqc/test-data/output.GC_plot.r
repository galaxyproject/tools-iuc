pdf("output.GC_plot.pdf")
gc=rep(c(60.78,41.18,47.06,56.86,29.41,27.45,37.25,78.43,58.82,50.98,49.02,62.75,68.63,54.90,52.94,35.29,43.14,39.22),times=c(3,3,5,7,1,2,2,1,1,3,2,1,1,1,3,1,2,1))
hist(gc,probability=T,breaks=100,xlab="GC content (%)",ylab="Density of Reads",border="blue",main="")
dev.off()
