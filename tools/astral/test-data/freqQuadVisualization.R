#!/usr/bin/env Rscript
red='#d53e4f';orange='#1d91c0';blue='#41b6c4';colormap = c(red,orange,blue)
require(reshape2);require(ggplot2);
dirPath = '.'; filePath = paste(dirPath,'/freqQuadCorrected.csv',sep=''); md<-read.csv(filePath,header=F,sep='\t'); md$value = md$V5/md$V6;
a<-length(levels(as.factor(md$V7)))*3.7; b<-4; sizes <- c(a,b);
md$V8<-reorder(md$V8,-md$value)
ggplot(data=md)+aes(x=V8,y=value,fill=V9)+geom_bar(stat='identity',color=1,width=0.8,position='dodge')+theme_bw()+theme(axis.text.x=element_text(angle=90))+scale_fill_manual(values=colormap,name='Topology')+geom_hline(yintercept=1/3,size=0.4,linetype=2)+ylab('relative freq.')+facet_wrap(~V7,scales='free_x')+xlab('')
pdfFile = paste(dirPath,'/relativeFreq.pdf',sep=''); ggsave(pdfFile,width = sizes[1], height= sizes[2]);
