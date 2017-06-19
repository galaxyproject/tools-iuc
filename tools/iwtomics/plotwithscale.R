if (require("IWTomics",character.only = TRUE,quietly = FALSE)) {
  args=commandArgs(TRUE)
  
  # get args names and values
  args_values=strsplit(args,'=')
  args_names=unlist(lapply(args_values,function(arg) arg[1]))
  names(args_values)=args_names
  args_values=lapply(args_values,function(arg) arg[2])
  # read filenames
  adjustedpvalue=args_values$adjustedpvalue
  iwtomicsrespdf=args_values$iwtomicsrespdf
  iwtomicssumpdf=args_values$iwtomicssumpdf
  iwtomicsrdata=args_values$iwtomicsrdata
  iwtomicstests=args_values$iwtomicstests
  iwtomicsselectedfeatures=args_values$iwtomicsselectedfeatures
  test_subset=paste0('c(',strsplit(args_values$test_subset,'\\|')[[1]],')')
  feature_subset=paste0('c(',strsplit(args_values$feature_subset,'\\|')[[1]],')')
  # read parameters (from test_subset on)
  i_scale_subset=which(args_names=='scale_subset')
  for(i in i_scale_subset:length(args)){
    eval(parse(text=args[[i]]))
  }
  
  # load RData
  load(iwtomicsrdata)
  # read testids and featureids and check them
  unlisted=lapply(seq_along(test_subset),
                   function(i){
                     test_subset_i=eval(parse(text=test_subset[i]))
                     feature_subset_i=eval(parse(text=feature_subset[i]))
                     test_subset_i=rep(test_subset_i,each=length(feature_subset_i))
                     feature_subset_i=rep(feature_subset_i,length.out=length(test_subset_i))
                     scale_subset_i=rep(scale_subset[i],length(test_subset_i))
                     return(list(test_subset=test_subset_i,feature_subset=feature_subset_i,scale_subset=scale_subset_i))
                   })
  test_subset=unlist(lapply(unlisted,function(l) l$test_subset))
  feature_subset=unlist(lapply(unlisted,function(l) l$feature_subset))
  scale_subset=unlist(lapply(unlisted,function(l) l$scale_subset))
  testids=as.character(read.delim(iwtomicstests,header=FALSE,sep='\t',stringsAsFactors=FALSE))
  featureids=as.character(read.delim(iwtomicsselectedfeatures,header=FALSE,sep='\t',stringsAsFactors=FALSE))
  id_features_subset=featureids[feature_subset]
  if(sum(testids!=paste(testInput(regionsFeatures_test)$id_region1,'vs',testInput(regionsFeatures_test)$id_region2))){
    write("Wrong test ids.", stderr())
    quit(save="no", status=10)
  }
  if(sum(featureids!=idFeatures(regionsFeatures_test))){
    write("Wrong feature ids.", stderr())
    quit(save="no", status=20)
  }
  # retrieve test and features_subset ids
  id_features_subset=featureids[feature_subset]
  if(sum(duplicated(paste0(test_subset,id_features_subset)))){
    write("Two scale thresholds selected for the same test and feature.", stderr())
    quit(save="no", status=30)
  }
  # If scale_subset=0, do not change the threshold
  default=(scale_subset==0)
  scale_subset=scale_subset[!default]
  test_subset=test_subset[!default]
  id_features_subset=id_features_subset[!default]
  
  # get scale threshold
  scale_threshold=lapply(regionsFeatures_test@test$result,
                         function(result) unlist(lapply(result,function(feature) feature$max_scale)))
  for(i in seq_along(test_subset)){
    if(scale_threshold[[test_subset[i]]][id_features_subset[i]]<scale_subset[i]){
      write("Scale threshold too high.", stderr())
      quit(save="no", status=40)
    }
    scale_threshold[[test_subset[i]]][id_features_subset[i]]=scale_subset[i]
  }
  
  # create adjustedvalue output
  pval=adjusted_pval(regionsFeatures_test,scale_threshold=scale_threshold)
  for(test in seq_along(pval)){
    for(id_feature in idFeatures(regionsFeatures_test)){
      write(paste0('Test: ',testids[test],', on feature ',id_feature),
              file=adjustedpvalue,append=TRUE,sep='\t')
        pval_i=as.data.frame(t(pval[[test]][[id_feature]]))
        row.names(pval_i)=paste('Scale',scale_threshold[[test]][[id_feature]])
        write.table(pval_i,file=adjustedpvalue,append=TRUE,sep='\t',quote=FALSE,row.names=TRUE,col.names=FALSE)
        write('',file=adjustedpvalue,append=TRUE,sep='\t')
      }
    }
  
  
  # plot test results
  pdf(iwtomicsrespdf,width=5,height=7)
  if(plottype=='boxplot'){
    # fix repeated probs
    probs=sort(unique(probs))
  }else{
    probs=c(0.25,0.5,0.75)
  }
  plotTest(regionsFeatures_test,alpha=testalpha,type=plottype,probs=probs,average=average,size=size,scale_threshold=scale_threshold,ask=FALSE)
  dev.off()
  
  # plot summary results
  if(groupby!='none'){
    tryCatch({
      pdf(iwtomicssumpdf,width=15,height=10)
      plotSummary(regionsFeatures_test,alpha=summaryalpha,only_significant=only_significant,groupby=groupby,scale_threshold=scale_threshold,ask=FALSE,append=TRUE)
      dev.off()
    }, error = function(err) {
      if (grepl('selected features with different resolution',err$message)) {
        write("Group by 'test' but selected features with different resolution.", stderr())
        quit(save="no", status=50) #error: groupby 'test' but selected features with different resolution.
      }
      write("Summary plot error. Please try again.", stderr())
      quit(save="no", status=60) #error 
    })
  }
  
}else{
  write("Missing IWTomics package. Please be sure to have it installed before using this tool.", stderr())
  quit(save="no", status=255)
}