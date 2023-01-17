# Load the IsoformSwitchAnalyzeR library
library(IsoformSwitchAnalyzeR, quietly = TRUE, warn.conflicts = FALSE)
library(GetoptLong, quietly = TRUE, warn.conflicts = FALSE)

# setup R error handling to go to stderr
options(
  show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

################################################################################
### Input Processing
################################################################################


# Collect arguments from command line
args <- commandArgs(trailingOnly = TRUE)

GetoptLong(
  "modeSelector=s","",
  "parentDir=s","",
  "readLength=i","",
  "annotation=s","",
  "transcriptome=s","",
  "fixStringTieAnnotationProblem","",
  "countFiles=s","",
  "toolSource=s","",
  "rObject=s","",
  "IFcutoff=f","",
  "geneExpressionCutoff=f","",
  "isoformExpressionCutoff=f","",
  "alpha=f","",
  "dIFcutoff=f","",
  "onlySigIsoforms","",
  "filterForConsequences","",
  "acceptedGeneBiotype=s","",
  "removeSingleIsformGenes","",
  "keepIsoformInAllConditions","",
  "correctForConfoundingFactors","",
  "overwriteIFvalues","",
  "reduceToSwitchingGenes","",
  "reduceFurtherToGenesWithConsequencePotential","",
  "keepIsoformInAllConditions2","",
  "minORFlength=i","",
  "orfMethod=s","",
  "PTCDistance=i","",
  "removeShortAAseq","",
  "removeLongAAseq","",
  "alsoSplitFastaFile","",
  "removeORFwithStop","",
  "onlySwitchingGenes",""
  #"analysisMode",       "H", 2, "character",
  #"genesToPLot",        "I", 2, "integer",
  #"gene",               "J", 2, "character",
  #"sortByQvals",        "K", 2, "logical",
  #"countGenes",         "M", 2, "logical",
  #"asFractionTotal",    "N", 2, "logical",
  #"plotGenes",          "O", 2, "logical",
  #"simplifyLocation",   "P", 2, "logical",
  #"removeEmptyConsequences",  "Q", 2, "logical",
  #"analysisOppositeConsequence",  "R", 2, "logical",
  #"pathToCPATresultFile", "S", 2, "character",
  #"pathToCPC2resultFile", "T", 2, "character",
  #"pathToPFAMresultFile", "U", 2, "character",
  #"pathToNetSurfP2resultFile", "V", 2, "character",
  #"pathToSignalPresultFile", "W", 2, "character",
  #"pathToIUPred2AresultFile","X", 2, "character",
  #"codingCutoff",       "Y", 2, "numeric",
  #"removeNoncodingORFs","Z", 2, "logical",
  #"minSignalPeptideProbability",  "aA", 2, "numeric",
  #"smoothingWindowSize","aB", 2, "integer",
  #"probabilityCutoff",  "aC", 2, "numeric",
  #"minIdrSize",         "aD", 2, "integer",
  #"annotateBindingSites","aE",2, "logical",
  #"minIdrBindingSize",  "aF", 2, "integer",
  #"minIdrBindingOverlapFrac", "aG", 2, "numeric",
  #"ntCutoff",           "aS", 2, "integer",
  #"ntFracCutoff",       "aH", 2, "numeric",
  #"ntJCsimCutoff",      "aI", 2, "numeric",
  #"AaCutoff",           "aJ", 2, "integer",
  #"AaFracCutoff",       "aK", 2, "numeric",
  #"AaJCsimCutoff",      "aM", 2, "numeric",
  #"removeNonConseqSwitches", "aN",  2,  "logical",
  #"rescaleTranscripts", "aO", 2, "logical",
  #"reverseMinus",       "aP", 2, "logical",
  #"addErrorbars",       "aQ", 2, "logical"
)

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "modeSelector",       "a", 1, "character",
  "parentDir",          "b", 2, "character",
  "readLength",         "c", 2, "integer",
  "annotation",         "d", 2, "character",
  "transcriptome",      "e", 2, "character",
  "fixStringTieAnnotationProblem", "f", 2, "logical",
  "countFiles",         "g", 2, "character",
  "toolSource",         "h", 2, "character",
  "rObject",            "i", 2, "character",
  "IFcutoff",           "j", 2, "numeric",
  "geneExpressionCutoff", "k", 2, "numeric",
  "isoformExpressionCutoff",  "m", 2, "numeric",
  "alpha",              "n", 2, "numeric",
  "dIFcutoff",          "o", 2, "numeric",
  "onlySigIsoforms",    "p", 2, "logical",
  "filterForConsequences",  "q", 2, "logical",
  "acceptedGeneBiotype","r", 2, "character",
  "removeSingleIsformGenes","s", 2, "logical",
  "keepIsoformInAllConditions", "t", 2, "logical",
  "correctForConfoundingFactors","u",2, "logical",
  "overwriteIFvalues",  "v", 2, "logical",
  "reduceToSwitchingGenes", "w", 2, "logical",
  "reduceFurtherToGenesWithConsequencePotential", "x", 2, "logical",
  "keepIsoformInAllConditions2", "y", 2, "logical",
  "minORFlength",       "z", 2, "integer",
  "orfMethod",          "A", 2, "character",
  "PTCDistance",        "B", 2, "integer",
  "removeShortAAseq",   "C", 2, "logical",
  "removeLongAAseq",    "D", 2, "logical",
  "alsoSplitFastaFile", "E", 2, "logical",
  "removeORFwithStop",  "F", 2, "logical",
  "onlySwitchingGenes", "G", 2, "logical"
  #"analysisMode",       "H", 2, "character",
  #"genesToPLot",        "I", 2, "integer",
  #"gene",               "J", 2, "character",
  #"sortByQvals",        "K", 2, "logical",
  #"countGenes",         "M", 2, "logical",
  #"asFractionTotal",    "N", 2, "logical",
  #"plotGenes",          "O", 2, "logical",
  #"simplifyLocation",   "P", 2, "logical",
  #"removeEmptyConsequences",  "Q", 2, "logical",
  #"analysisOppositeConsequence",  "R", 2, "logical",
  #"pathToCPATresultFile", "S", 2, "character",
  #"pathToCPC2resultFile", "T", 2, "character",
  #"pathToPFAMresultFile", "U", 2, "character",
  #"pathToNetSurfP2resultFile", "V", 2, "character",
  #"pathToSignalPresultFile", "W", 2, "character",
  #"pathToIUPred2AresultFile","X", 2, "character",
  #"codingCutoff",       "Y", 2, "numeric",
  #"removeNoncodingORFs","Z", 2, "logical",
  #"minSignalPeptideProbability",  "aA", 2, "numeric",
  #"smoothingWindowSize","aB", 2, "integer",
  #"probabilityCutoff",  "aC", 2, "numeric",
  #"minIdrSize",         "aD", 2, "integer",
  #"annotateBindingSites","aE",2, "logical",
  #"minIdrBindingSize",  "aF", 2, "integer",
  #"minIdrBindingOverlapFrac", "aG", 2, "numeric",
  #"ntCutoff",           "aS", 2, "integer",
  #"ntFracCutoff",       "aH", 2, "numeric",
  #"ntJCsimCutoff",      "aI", 2, "numeric",
  #"AaCutoff",           "aJ", 2, "integer",
  #"AaFracCutoff",       "aK", 2, "numeric",
  #"AaJCsimCutoff",      "aM", 2, "numeric",
  #"removeNonConseqSwitches", "aN",  2,  "logical",
  #"rescaleTranscripts", "aO", 2, "logical",
  #"reverseMinus",       "aP", 2, "logical",
  #"addErrorbars",       "aQ", 2, "logical"
  ),
byrow = TRUE, ncol = 4
)


opt <- getopt(spec)

# Run IsoformSwitchAnalyzeR

if (opt$modeSelector == 'data_import') {

  # Analysis part one
  ###################

  quantification_data <- importIsoformExpression(
    parentDir = opt$parentDir,
    addIsofomIdAsColumn = TRUE, #not optional, required to be TRUE
    readLength = opt$readLength
  )

  ### Make design matrix
  myDesign <- data.frame(
    sampleID = colnames(quantification_data$abundance)[-1],
    condition = gsub('[[:digit:]]+', '', colnames(quantification_data$abundance)[-1])
  )

  if (opt$toolSource == 'stringtie'){
    SwitchList <- importRdata(
      isoformCountMatrix   = quantification_data$counts,
      isoformRepExpression = quantification_data$abundance,
      designMatrix         = myDesign,
      isoformExonAnnoation = opt$annotation,
      isoformNtFasta       = opt$transcriptome,
      showProgress = TRUE,
      fixStringTieAnnotationProblem = opt$fixStringTieAnnotationProblem)
  }else{
    SwitchList <- importRdata(
      isoformCountMatrix   = quantification_data$counts,
      isoformRepExpression = quantification_data$abundance,
      designMatrix         = myDesign,
      isoformExonAnnoation = opt$annotation,
      isoformNtFasta       = opt$transcriptome,
      showProgress = TRUE)
  }
  

  geneCountMatrix <- extractGeneExpression(
    SwitchList,
    extractCounts = TRUE,
    addGeneNames = FALSE,
    addIdsAsColumns = FALSE
  )

  if (opt$countFiles == 'deseq2') {

    # Generate count files for DESeq2

    l<-as.list(as.data.frame(geneCountMatrix))

    sample_names <- colnames(as.data.frame(geneCountMatrix))
    gene_names <- row.names(as.data.frame(geneCountMatrix))

    for (index in 1:length(l)){
      tabular_expression <- data.frame(gene_names,l[index])
      colnames(tabular_expression) <- c("Geneid",sample_names[index])
      filename <- paste(sample_names[index],"dataset.tabular", sep="_")
      output_path <- paste("./count_files/",filename,sep="")
      write.table(tabular_expression,output_path, sep="\t",row.names=FALSE,quote = FALSE)
    }
    }else if (opt$countFiles == 'edger') {
      expressionDF <- data.frame(geneCountMatrix)
      gene_names <- row.names(expressionDF)
      expressionDF <- cbind(gene_names,expressionDF)
      write.table(as.data.frame(expressionDF),"./count_files/edgeR.tabular", sep="\t",row.names=FALSE,quote = FALSE)
    } 

  save(SwitchList,file="SwitchList.Rda")

} 

if (opt$modeSelector == 'first_step') {

  # Second part of the analysis
  #############################

  load(file=opt$rObject)

  ### Filter
  SwitchList <- preFilter( SwitchList,
                          IFcutoff=opt$IFcutoff,
                          geneExpressionCutoff = opt$geneExpressionCutoff,
                          isoformExpressionCutoff = opt$isoformExpressionCutoff,
                          removeSingleIsoformGenes = opt$removeSingleIsformGenes,
                          #reduceToSwitchingGenes=FALSE,
                          onlySigIsoforms = opt$onlySigIsoforms,
                          keepIsoformInAllConditions=opt$keepIsoformInAllConditions,
                          alpha=opt$alpha,
                          dIFcutoff = opt$dIFcutoff,
                          acceptedGeneBiotype = opt$acceptedGeneBiotype
  )

  ### Test for isoform switches
  SwitchList <- isoformSwitchTestDEXSeq( SwitchList,
                                        alpha = opt$alpha,
                                        dIFcutoff = opt$dIFcutoff,
                                        ### Advanced arguments
                                        correctForConfoundingFactors=opt$correctForConfoundingFactors,
                                        overwriteIFvalues=opt$overwriteIFvalues,
                                        reduceToSwitchingGenes = opt$reduceToSwitchingGenes,
                                        reduceFurtherToGenesWithConsequencePotential = opt$reduceFurtherToGenesWithConsequencePotential,
                                        onlySigIsoforms = opt$onlySigIsoforms,
                                        keepIsoformInAllConditions = opt$keepIsoformInAllConditions2,
                                        showProgress = TRUE,
  ) 

  ### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
  #SwitchList <- addORFfromGTF( SwitchList,
  #                             pathToGTF = ,
  #                             ### Advanced argument
  #                             overwriteExistingORF = FALSE,
  #                             onlyConsiderFullORF = FALSE,
  #                             removeNonConvensionalChr = FALSE,
  #                             PTCDistance = 50,)

  SwitchList <- analyzeNovelIsoformORF( SwitchList,
                                        analysisAllIsoformsWithoutORF = TRUE, # also analyse all those annoatated as without CDS in ref annottaion
                                        ### Advanced argument
                                        minORFlength = opt$minORFlength,
                                        orfMethod = opt$orfMethod,
                                        PTCDistance = opt$PTCDistance,
                                        startCodons = "ATG",
                                        stopCodons = c("TAA", "TAG", "TGA"),
                                        showProgress = TRUE,
  )

  ### Extract Sequences
  SwitchList <- extractSequence( SwitchList,
                                onlySwitchingGenes = opt$onlySwitchingGenes,
                                alpha = opt$alpha,
                                dIFcutoff = opt$dIFcutoff,
                                extractNTseq = TRUE,
                                extractAAseq = TRUE,
                                removeShortAAseq = opt$removeShortAAseq,
                                removeLongAAseq = opt$removeLongAAseq,
                                alsoSplitFastaFile = opt$alsoSplitFastaFile,
                                removeORFwithStop=opt$removeORFwithStop,
                                addToSwitchAnalyzeRlist = TRUE,
                                writeToFile = TRUE,
                                pathToOutput = getwd(),
                                outputPrefix='isoformSwitchAnalyzeR_isoform',
                                forceReExtraction = FALSE,
                                quiet=FALSE
  )

  ### Summary

  switchSummary <- extractSwitchSummary( SwitchList,
                        filterForConsequences=opt$filterForConsequences,
                        alpha=opt$alpha,
                        dIFcutoff = opt$dIFcutoff,
                        onlySigIsoforms = opt$onlySigIsoforms,
                        # includeCombined=nrow(unique(switchAnalyzeRlist$isoformFeatures[,c('condition_1','condition_1')]))
  )

  save(SwitchList,file="SwitchList.Rda")
  write.table(switchSummary, file="switchSummary.tsv",quote=FALSE, sep='\t', col.names = TRUE, row.names=FALSE)

}

if (opt$modeSelector == 'second_step'){
  print("Hello, God!")


  # Second part of the analysis
  #############################

  load(file=opt$rObject)


  ### Add annotation

  if (!is.null(opt$pathToCPATresultFile)) {
    SwitchList <- analyzeCPAT( SwitchList,
                                pathToCPATresultFile = opt$pathToCPATresultFile,
                                codingCutoff = opt$codingCutoff,
                                removeNoncodinORFs        = opt$removeNoncodingORFs
    )
  }

  if (!is.null(opt$pathToCPC2resultFile)) {
    SwitchList <- analyzeCPC2( SwitchList,
                              pathToCPC2resultFile = opt$pathToCPC2resultFile,
                              removeNoncodinORFs = opt$removeNoncodingORFs # because ORF was predicted de Novo
    )
  }

  if (!is.null(opt$pathToPFAMresultFile)) {

    pfamFiles <- list.files(path=opt$pathToPFAMresultFile,
                        full.names=TRUE
                        )

    SwitchList <- analyzePFAM( SwitchList,
                              pathToPFAMresultFile =  pfamFiles,
                              showProgress=FALSE
    )
  }

  if (!is.null(opt$pathToNetSurfP2resultFile)) {

    netsurfFiles <- list.files(path=opt$pathToNetSurfP2resultFile,
                    full.names=TRUE
    )

    SwitchList <- analyzeNetSurfP2( SwitchList,
                                    pathToNetSurfP2resultFile =  netsurfFiles,
                                    smoothingWindowSize = opt$smoothingWindowSize,
                                    probabilityCutoff = opt$probabilityCutoff,
                                    minIdrSize = opt$minIdrSize,
                                    showProgress = TRUE
    )
  }

  if (!is.null(opt$pathToIUPred2AresultFile)) {

    SwitchList <- analyzeIUPred2A(SwitchList,
                                  pathToIUPred2AresultFile = opt$pathToIUPred2AresultFile,
                                  smoothingWindowSize = opt$smoothingWindowSize,
                                  probabilityCutoff = opt$probabilityCutoff,
                                  minIdrSize = opt$minIdrSize,
                                  annotateBindingSites = opt$annotateBindingSites,
                                  minIdrBindingSize = opt$minIdrSize,
                                  minIdrBindingOverlapFrac = opt$minIdrBindingOverlapFrac,
                                  showProgress = TRUE,
                                  quiet = FALSE
                                  )
  }

  if (!is.null(opt$pathToSignalPresultFile)) {

    signalpFiles <- list.files(path=opt$pathToSignalPresultFile,
                            full.names=TRUE
    )

    SwitchList <- analyzeSignalP( SwitchList,
                                  pathToSignalPresultFile = signalpFiles,
                                  minSignalPeptideProbability = opt$minSignalPeptideProbability
    )
  }

  SwitchList <- analyzeAlternativeSplicing( SwitchList,
                                            onlySwitchingGenes=opt$onlySwitchingGenes,
                                            alpha=opt$alpha,
                                            dIFcutoff = opt$dIFcutoff,
                                            showProgress=TRUE
  )



  SwitchList <- analyzeIntronRetention( SwitchList,
                                        onlySwitchingGenes = opt$onlySwitchingGenes,
                                        alpha = opt$alpha,
                                        dIFcutoff = opt$dIFcutoff,
                                        showProgress = TRUE
  )
  
  SwitchList <- analyzeSwitchConsequences( SwitchList,
                                          consequencesToAnalyze='all',
                                            # c(
                                            # 'intron_retention',
                                            # 'coding_potential',
                                            # 'ORF_seq_similarity',
                                            # 'NMD_status',
                                            # 'domains_identified',
                                            #'IDR_identified',
                                            #'IDR_type',
                                            # 'signal_peptide_identified'),
                                          alpha=opt$alpha,
                                          dIFcutoff=opt$dIFcutoff,
                                          onlySigIsoforms=opt$onlySigIsoforms,
                                          ntCutoff=opt$ntCutoff,
                                          ntJCsimCutoff=0.8,
                                          AaCutoff=10,
                                          AaFracCutoff=0.5,
                                          AaJCsimCutoff=0.9,
                                          removeNonConseqSwitches=TRUE,
                                          showProgress=TRUE,
                                          if (!missing(ntFracCutoff)){
                                            ntFracCutoff=opt$ntFracCutoff
                                          }

  )

  ### Visual analysis
  # Top genes

  if (analysisMode == 'single'){
      ## This function enables a full analysis of a specific gene containing an isoform switch
      output_file <- file.path(getwd(), "single_gene.pdf")
      pdf(file = output_file, onefile = FALSE, height=6, width = 9)
      switchPlot(
        ### Core arguments
        SwitchList,
        gene = opt$gene,
        #isoform_id = NULL,
        condition1 = "cancer", # NEEDS TO BE CHANGED!!!!
        condition2 = "health", # NEEDS TO BE CHANGED!!!!
        ### Advanced arguments
        IFcutoff = opt$IFcutoff,
        dIFcutoff = opt$dIFcutoff,
        # alphas = c(0.05, 0.001),
        rescaleTranscripts = opt$rescaleTranscripts,
        reverseMinus = opt$rescaleTranscripts,
        addErrorbars = opt$addErrorbars,
        logYaxis = opt$logYaxis,
        localTheme = theme_bw(base_size = 8),
      )
      dev.off()

  }else{

    mostSwitchingGene <- extractTopSwitches( ## options included in plot_mode
      SwitchList,
      n = Inf,
      filterForConsequences=opt$filterForConsequences,
      extractGenes=TRUE,
      alpha=opt$alpha,
      dIFcutoff = opt$dIFcutoff,
      inEachComparison=FALSE, # not included
      sortByQvals=opt$sortByQvals
    )

    switchPlotTopSwitches( SwitchList,
                          alpha = opt$alpha,
                          dIFcutoff = opt$dIFcutoff,
                          onlySigIsoforms = opt$onlySigIsoforms,
                          n= opt$genesToPLot
                          sortByQvals=opt$sortByQvals,
                          filterForConsequences = opt$filterForConsequences,
                          pathToOutput = getwd(),
                          splitComparison=opt$splitComparison,
                          splitFunctionalConsequences = opt$splitFunctionalConsequences,
                          IFcutoff=opt$IFcutoff,
                          fileType = "png",
    )
    
    output_file <- file.path(getwd(), "extractConsequencesSummary.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 10)
    consequenceSummary <- extractConsequenceSummary(
      SwitchList,
      consequencesToAnalyze='all',
      includeCombined=FALSE, # not included
      asFractionTotal=opt$asFractionTotal,
      alpha=opt$asFractionTotal,
      dIFcutoff=0.1,
      plot=TRUE,
      plotGenes=FALSE,
      simplifyLocation = TRUE,
      removeEmptyConsequences = FALSE,
      returnResult=TRUE, # not included
      localTheme=theme_bw()
    )
    dev.off()

    output_file <- file.path(getwd(), "consequencesEnrichment.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    consequenceEnrichment <- extractConsequenceEnrichment(
      SwitchList,
      consequencesToAnalyze = 'all',
      alpha=0.05,
      dIFcutoff = 0.1,
      countGenes = TRUE,
      analysisOppositeConsequence=FALSE,
      plot=TRUE, #not included
      localTheme = theme_bw(base_size = 12), #not included
      minEventsForPlotting = 10, #not included
      returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
    )
    dev.off()

    output_file <- file.path(getwd(), "splicingEnrichment.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    sprincingEnrichment <- extractSplicingEnrichment(
      SwitchList,
      splicingToAnalyze = 'all',
      alpha = 0.05,
      dIFcutoff = 0.1,
      onlySigIsoforms = FALSE,
      countGenes = TRUE,
      plot = TRUE,
      minEventsForPlotting = 10, #not included
      returnResult = TRUE # not included
    )
    dev.off()

    output_file <- file.path(getwd(), "splicingSummary.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    splicingSummary <- extractSplicingSummary(
      SwitchList,
      splicingToAnalyze = 'all',
      asFractionTotal = opt$asFractionTotal,
      alpha = opt$alpha,
      dIFcutoff = opt$dIFcutoff,
      onlySigIsoforms = opt$onlySigIsoforms,
      plot = TRUE,
      plotGenes = FALSE,
      localTheme = theme_bw(),
      returnResult = TRUE
    )
    dev.off()

    ### Volcano like plot:
    output_file <- file.path(getwd(), "volcanoPlot.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    ggplot(data=SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
      geom_point(
        aes( color=abs(dIF) > opt$dIFcutoff & isoform_switch_q_value < opt$alpha ), # default cutoff
        size=1
      ) +
      geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
      geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
      facet_wrap( ~ condition_2) +
      #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
      scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
      labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
      theme_bw()
    dev.off()


    ### Switch vs Gene changes:
    output_file <- file.path(getwd(), "switchGene.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    ggplot(data=SwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
      geom_point(
        aes( color=abs(dIF) > opt$dIFcutoff & isoform_switch_q_value < opt$alpha ), # default cutoff
        size=1
      ) + 
      facet_wrap(~ condition_2) +
      #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
      geom_hline(yintercept = 0, linetype='dashed') +
      geom_vline(xintercept = 0, linetype='dashed') +
      scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
      labs(x='Gene log2 fold change', y='dIF') +
      theme_bw()
    dev.off()

    output_file <- file.path(getwd(), "splicingGenomewide.pdf")
    pdf(file = output_file, onefile = FALSE, height=6, width = 9)
    splicingGenomeWide <- extractSplicingGenomeWide(
      SwitchList,
      featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
      splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
      plot=TRUE,
      returnResult=TRUE  # Preventing the summary statistics to be returned as a data.frame
    )
    dev.off()
  }
}