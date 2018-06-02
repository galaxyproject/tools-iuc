if (require("IWTomics",character.only = TRUE,quietly = FALSE)) {
  args=commandArgs(TRUE)
  
  # get args names and values
  args_values=strsplit(args,'=')
  args_names=unlist(lapply(args_values,function(arg) arg[1]))
  names(args_values)=args_names
  args_values=lapply(args_values,function(arg) arg[2])
  # read filenames
  adjustedpvaluematrix=args_values$adjustedpvaluematrix
  iwtomicsrespdf=args_values$iwtomicsrespdf
  iwtomicssumpdf=args_values$iwtomicssumpdf
  regionids=args_values$regionids
  featureids=args_values$featureids
  rdatafile=args_values$rdatafile
  iwtomicsrdata=args_values$iwtomicsrdata
  iwtomicstests=args_values$iwtomicstests
  iwtomicsselectedfeatures=args_values$iwtomicsselectedfeatures
  # read parameters (from region1 on)
  i_region1=which(args_names=='region1')
  for(i in i_region1:length(args)){
    eval(parse(text=args[[i]]))
  }
  
  # load RData
  load(rdatafile)
  # read regionids and featureids
  regionids=as.character(read.delim(regionids,header=FALSE,sep='\t',stringsAsFactors=FALSE))
  featureids=as.character(read.delim(featureids,header=FALSE,sep='\t',stringsAsFactors=FALSE))
  # retrieve region1, region2 and features_subset ids and check they are in the RData
  id_region1=regionids[region1]
  id_region2=regionids[region2]
  id_features_subset=featureids[features_subset]
  if(length(setdiff(c(id_region1,id_region2),idRegions(regionsFeatures)))!=0){
    write("Wrong region ids.", stderr())
    quit(save="no", status=10)
  }
  if(length(setdiff(id_features_subset,idFeatures(regionsFeatures)))!=0){
    write("Wrong feature ids.", stderr())
    quit(save="no", status=20)
  }
  if(sum(duplicated(paste0(id_region1,id_region2)))){
    write("Same test repeated multiple times.", stderr())
    quit(save="no", status=30)
  }
  
  # perform test
  tryCatch({
    # fix repeated probs
    if(statistics=='quantile'){
      # fix repeated probs
      testprobs=sort(unique(testprobs))
    }else{
      testprobs=0.5
    }
    regionsFeatures_test=IWTomicsTest(regionsFeatures,id_region1,id_region2,id_features_subset,
                                      statistics=statistics,probs=testprobs,B=B)
    # create adjustedvaluematrix output
    for(test in seq_along(id_region1)){
      for(id_feature in id_features_subset){
        write(paste0('Test: ',id_region1[test],' vs ',id_region2[test],', on feature ',id_feature),
              file=adjustedpvaluematrix,append=TRUE,sep='\t')
        pval=regionsFeatures_test@test$result[[test]][[id_feature]]$adjusted_pval_matrix
        row.names(pval)=paste('Scale',rev(seq_len(nrow(pval))))
        write.table(pval,file=adjustedpvaluematrix,append=TRUE,sep='\t',quote=FALSE,row.names=TRUE,col.names=FALSE)
        write('',file=adjustedpvaluematrix,append=TRUE,sep='\t')
      }
    }
  }, error = function(err) {
    write("Testing error.", stderr())
    quit(save="no", status=40) #error testing
  })
  
  # plot test results
  pdf(iwtomicsrespdf,width=5,height=7)
  if(plottype=='boxplot'){
    # fix repeated probs
    probs=sort(unique(probs))
  }else{
    probs=c(0.25,0.5,0.75)
  }
  plotTest(regionsFeatures_test,alpha=testalpha,type=plottype,probs=probs,average=average,size=size,ask=FALSE)
  dev.off()
  
  # plot summary results
  if(groupby!='none'){
    pdf(iwtomicssumpdf,width=15,height=10)
    plotSummary(regionsFeatures_test,alpha=summaryalpha,only_significant=only_significant,groupby=groupby,ask=FALSE,append=TRUE)
    dev.off()
  }
  
  # create output
  write.table(as.data.frame(t(paste(id_region1,'vs',id_region2))),file=iwtomicstests,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  write.table(as.data.frame(t(idFeatures(regionsFeatures_test))),file=iwtomicsselectedfeatures,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  save(regionsFeatures_test,file=iwtomicsrdata)
}else{
  write("Missing IWTomics package. Please be sure to have it installed before using this tool.", stderr())
  quit(save="no", status=255)
}