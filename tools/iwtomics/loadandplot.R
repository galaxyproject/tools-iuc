if (require("IWTomics",character.only = TRUE,quietly = FALSE)) {
  args=commandArgs(TRUE)
  
  # get args names and values
  args_values=strsplit(args,'=')
  args_names=unlist(lapply(args_values,function(arg) arg[1]))
  names(args_values)=args_names
  args_values=lapply(args_values,function(arg) arg[2])
  # read filenames
  outrdata=args_values$outrdata
  outregions=args_values$outregions
  outfeatures=args_values$outfeatures
  outpdf=args_values$outpdf
  regionspaths=unlist(strsplit(args_values$regionspaths,'\\|'))
  if("regionsheaderfile" %in% args_names){
    # the file regionsheaderfile must contain as first column the (unique) regionsfilenames, 
    # as second column the corresponding ids and as third column the names
    tryCatch({
      regionsheader=read.delim(args_values$regionsheaderfile,header=FALSE,stringsAsFactors=FALSE,row.names=1,sep="\t")
      regionsfilenames=unlist(strsplit(args_values$regionsfilenames,'\\|'))
      if(length(setdiff(regionsfilenames,row.names(regionsheader)))) {
        write("IWTomics message: Not all region files are present in the first column of header file for regions.", stderr())
        quit(save="no", status=11)
      }
      id_regions=regionsheader[regionsfilenames,1]
      name_regions=regionsheader[regionsfilenames,2]
    }, error = function(err) {
      write("IWTomics message: An error has occurred reading the header file for regions. Please try again.", stderr())
      quit(save="no", status=10) #error on header file
    })
  }else{
    eval(parse(text=args[[which(args_names=='regionsgalaxyids')]]))
    id_regions=paste0('data_',regionsgalaxyids)
    name_regions=paste0('data_',regionsgalaxyids)
  }
  featurespaths=unlist(strsplit(args_values$featurespaths,'\\|'))
  if("featuresheaderfile" %in% args_names){
    # the file featuresheaderfile must contain as first column the (unique) featuresfilenames, 
    # as second column the corresponding ids and as third column the names
    tryCatch({
      featuresheader=read.delim(args_values$featuresheaderfile,header=FALSE,stringsAsFactors=FALSE,row.names=1,sep="\t")
      featuresfilenames=unlist(strsplit(args_values$featuresfilenames,'\\|'))
      if(length(setdiff(featuresfilenames,row.names(featuresheader)))) {
        write("IWTomics message: Not all feature files are present in the first column of header file for features.", stderr())
        quit(save="no", status=21)
      }
      id_features=featuresheader[featuresfilenames,1]
      name_features=featuresheader[featuresfilenames,2]
    }, error = function(err) {
      write("IWTomics message: An error has occurred reading the header file for features. Please try again.", stderr())
      quit(save="no", status=20) #error on header file
    })
  }else{
    eval(parse(text=args[[which(args_names=='featuresgalaxyids')]]))
    id_features=paste0('data_',featuresgalaxyids)
    name_features=paste0('data_',featuresgalaxyids)
  }
  # read parameters (from smoothing on)
  i_smoothing=which(args_names=='smoothing')
  for(i in i_smoothing:length(args)){
    eval(parse(text=args[[i]]))
  }
  
  # load data
  tryCatch({
    regionsFeatures=IWTomicsData(regionspaths,featurespaths,alignment,
                                id_regions,name_regions,id_features,name_features,start.are.0based=start.are.0based)
  }, error = function(err) {
    if(grepl('invalid format',err$message)){
      write("IWTomics message: Not enough columns in input file.", stderr())
      quit(save="no", status=31) # error, not enough columns in input file
    }else if(grepl('duplicated regions',err$message)){
      write("IWTomics message: Duplicated regions in region file.", stderr())
      quit(save="no", status=32) # error, duplicated regions in region file
    }else if(grepl('duplicated windows',err$message)){
      write("IWTomics message: Duplicated windows in feature file.", stderr())
      quit(save="no", status=33) # error, duplicated windows in feature file
    }else if(grepl('overlapping windows',err$message)){
      write("IWTomics message: Overlapping windows in feature file.", stderr())
      quit(save="no", status=34) # error, overlapping windows in feature file
    }else if(grepl('not all regions in datasets',err$message)){
      write("IWTomics message: Windows in feature files do not cover all regions in region files.", stderr())
      quit(save="no", status=35) # error, windows in feature files do not cover all regions in region files
    }else if(grepl('ifferent size windows',err$message)){
      write("IWTomics message: All windows in a feature file must have the same size.", stderr())
      quit(save="no", status=36) # error, all windows in a feature files must have the same size
    }
    #error loading data
    write("IWTomics message: An error has occurred reading the data. Please try again.", stderr())
    quit(save="no", status=30)    
 })
  
  # smooth data
  if(smoothing!='no'){
    tryCatch({
      if(smoothing=='locpoly'){
        dist_knots=10
      }else if(smoothing=='kernel'){
        degree=3
        dist_knots=10
      }else if(smoothing=='splines'){
        bandwidth=5
      }
      if(alignment=='scale'){
        if(scale==0){
          regionsFeatures=smooth(regionsFeatures,type=smoothing,fill_gaps=fill_gaps,
                                bandwidth=bandwidth,degree=degree,dist_knots=dist_knots)
        }else{
          regionsFeatures=smooth(regionsFeatures,type=smoothing,fill_gaps=fill_gaps,
                                bandwidth=bandwidth,degree=degree,dist_knots=dist_knots,scale_grid=scale)
        }
      }else{
        regionsFeatures=smooth(regionsFeatures,type=smoothing,fill_gaps=fill_gaps,
                              bandwidth=bandwidth,degree=degree,dist_knots=dist_knots)
      }
    }, error = function(err) {
      write("IWTomics message: An error has occurred smoothing the data. Please try again.", stderr())
      quit(save="no", status=40) #error on smoothing
    })
  }
  
  # plot data
  pdf(outpdf,width=10,height=8)
  if(plottype=='boxplot'){
    # fix repeated probs
    probs=sort(unique(probs))
  }else{
    probs=c(0.25,0.5,0.75)
  }
  plot(regionsFeatures,type=plottype,probs=probs,average=average,size=size,ask=FALSE)
  dev.off()
  
  # create output
  #write.table(cbind(unlist(strsplit(args_values$regionsfilenames,'\\|')),idRegions(regionsFeatures),nameRegions(regionsFeatures)),
              #file=outregions,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  write.table(as.data.frame(t(idRegions(regionsFeatures))),file=outregions,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  #write.table(cbind(unlist(strsplit(args_values$featuresfilenames,'\\|')),idFeatures(regionsFeatures),nameFeatures(regionsFeatures)),
              #file=outfeatures,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  write.table(as.data.frame(t(idFeatures(regionsFeatures))),file=outfeatures,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
  save(regionsFeatures,file=outrdata)
}else{
  write("IWTomics message: Missing IWTomics package. Please be sure to have it installed before using this tool.", stderr())
  quit(save="no", status=255)
}