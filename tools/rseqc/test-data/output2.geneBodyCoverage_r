pairend_strandspecific_51mer_hg19_chr1_1_100000_bam <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.1 <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.2 <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
data_matrix <- matrix(c(pairend_strandspecific_51mer_hg19_chr1_1_100000_bam,pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.1,pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.2), byrow=T, ncol=100)
rowLabel <- c("pairend_strandspecific_51mer_hg19_chr1_1_100000_bam","pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.1","pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.2")


pdf("output.geneBodyCoverage.heatMap.pdf")
rc <- cm.colors(ncol(data_matrix))
heatmap(data_matrix, scale=c("none"),keep.dendro=F, labRow = rowLabel ,Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),margins = c(6, 8),ColSideColors = rc,cexRow=1,cexCol=1,xlab="Gene body percentile (5'->3')", add.expr=x_axis_expr <- axis(side=1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("1","10","20","30","40","50","60","70","80","90","100")))
dev.off()


pdf("output.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(3)
plot(x,pairend_strandspecific_51mer_hg19_chr1_1_100000_bam,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
lines(x,pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.1,type='l',col=icolor[2])
lines(x,pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.2,type='l',col=icolor[3])
legend(0,1,fill=icolor[1:3], legend=c('pairend_strandspecific_51mer_hg19_chr1_1_100000_bam','pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.1','pairend_strandspecific_51mer_hg19_chr1_1_100000_bam.2'))
dev.off()
