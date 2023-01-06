# Load the IsoformSwitchAnalyzeR library
library(IsoformSwitchAnalyzeR, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)

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

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "modeSelector",   "a",  1,  "character",
  "parentDir",      "b",  2,  "character",
  "readLength",     "c",  2,  "integer",
  "annotation",     "d",  2,  "character",
  "transcriptome",  "e",  2,  "character",
  "fixStringTieAnnotationProblem",  "f",  2,  "logical",
  "countFiles",     "g",  2,  "character"),
byrow = TRUE, ncol = 4
)
opt <- getopt(spec)

# Run IsoformSwitchAnalyzeR

if (opt$modeSelector == 'data_import') {

  # Analysis part one
  ###################

  quantification_data <- importIsoformExpression(
    parentDir = opt$parentDir,
    addIsofomIdAsColumn = TRUE,
    readLength = opt$readLength
  )

  ### Make design matrix
  myDesign <- data.frame(
    sampleID = colnames(quantification_data$abundance)[-1],
    condition = gsub('[[:digit:]]+', '', colnames(quantification_data$abundance)[-1])
  )

  SwitchList <- importRdata(
    isoformCountMatrix   = quantification_data$counts,
    isoformRepExpression = quantification_data$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = opt$annotation,
    isoformNtFasta       = opt$transcriptome,
    fixStringTieAnnotationProblem = opt$fixStringTieAnnotationProblem,
    showProgress = TRUE
  )

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

  ### Filter
  SwitchList <- preFilter( SwitchList,
                          IFcutoff=0.01,
                          geneExpressionCutoff = 1,
                          isoformExpressionCutoff = 0,
                          #acceptedGeneBiotype = NULL,
                          removeSingleIsoformGenes = TRUE,
                          #reduceToSwitchingGenes=FALSE,
                          onlySigIsoforms = FALSE,
                          keepIsoformInAllConditions=FALSE,
                          alpha=0.05,
                          dIFcutoff = 0.1,
  )

  ### Test for isoform switches
  SwitchList <- isoformSwitchTestDEXSeq( SwitchList,
                                        alpha = 0.05,
                                        dIFcutoff = 0.1,
                                        ### Advanced arguments
                                        correctForConfoundingFactors=TRUE,
                                        overwriteIFvalues=TRUE,
                                        reduceToSwitchingGenes = TRUE,
                                        reduceFurtherToGenesWithConsequencePotential = TRUE,
                                        onlySigIsoforms = FALSE,
                                        keepIsoformInAllConditions = TRUE,
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
                                        #genomeObject = NULL,
                                        ### Advanced argument
                                        minORFlength = 100,
                                        orfMethod = 'longest.AnnotatedWhenPossible',
                                        PTCDistance = 50,
                                        startCodons = "ATG",
                                        stopCodons = c("TAA", "TAG", "TGA"),
                                        showProgress = TRUE,
  )

  ### Extract Sequences
  SwitchList <- extractSequence( SwitchList,
                                onlySwitchingGenes = TRUE,
                                alpha = 0.05,
                                dIFcutoff = 0.1,
                                extractNTseq = TRUE,
                                extractAAseq = TRUE,
                                removeShortAAseq = TRUE,
                                removeLongAAseq = FALSE,
                                alsoSplitFastaFile = FALSE,
                                removeORFwithStop=TRUE,
                                addToSwitchAnalyzeRlist = TRUE,
                                writeToFile = TRUE,
                                pathToOutput = getwd(),
                                outputPrefix='isoformSwitchAnalyzeR_isoform',
                                forceReExtraction = FALSE,
                                quiet=FALSE
  )

  ### Summary

  switchSummary <- extractSwitchSummary( SwitchList,
                        filterForConsequences=FALSE,
                        alpha=0.05,
                        dIFcutoff = 0.1,
                        onlySigIsoforms = FALSE,
                        # includeCombined=nrow(unique(switchAnalyzeRlist$isoformFeatures[,c('condition_1','condition_1')]))
  )

  switchSummary
}

  if (opt$modeSelector == 'second_step') {

  # Second part of the analysis
  #############################

  ### Add annotation

  SwitchList <- analyzeCPAT( SwitchList,
                              pathToCPATresultFile = "/home/laptop/stringtie_test/isoformSwitch_outputs/result_cpat.txt",
                              codingCutoff = 0.725,
                              removeNoncodinORFs        = TRUE
  )

  SwitchList <- analyzeCPC2( SwitchList,
                            pathToCPC2resultFile = "/home/laptop/stringtie_test/isoformSwitch_outputs/CPC2_standalone-1.0.1/result_cpc2.txt",
                            removeNoncodinORFs = TRUE # because ORF was predicted de Novo
  )


  SwitchList <- analyzePFAM( SwitchList,
                            pathToPFAMresultFile =  "/home/laptop/stringtie_test/isoformSwitch_outputs/pfam_results.txt",
                            showProgress=FALSE
  )

  SwitchList <- analyzeNetSurfP2( SwitchList,
                                  pathToNetSurfP2resultFile =  "/home/laptop/stringtie_test/isoformSwitch_outputs/networksurfp2_results.txt",
                                  smoothingWindowSize = 5,
                                  probabilityCutoff = 0.5,
                                  minIdrSize = 30,
                                  showProgress = TRUE,
  )

  SwitchList <- analyzeAlternativeSplicing( SwitchList,
                                            onlySwitchingGenes=FALSE,
                                            alpha=0.05,
                                            dIFcutoff = 0.1,
                                            showProgress=TRUE
  )



  SwitchList <- analyzeSignalP( SwitchList,
                                pathToSignalPresultFile = c("/home/laptop/stringtie_test/isoformSwitch_outputs/signalp_files/first_group/output_protein_type.txt",
                                                            "/home/laptop/stringtie_test/isoformSwitch_outputs/signalp_files/second_group/output_protein_type.txt"),
                                minSignalPeptideProbability = 0.5,
  )

  SwitchList <- analyzeIntronRetention( SwitchList,
                                        onlySwitchingGenes = TRUE,
                                        alpha = 0.05,
                                        dIFcutoff = 0.1,
                                        showProgress = TRUE
  )

  SwitchList <- analyzeSwitchConsequences( SwitchList,
                                          consequencesToAnalyze=c(
                                            'intron_retention',
                                            'coding_potential',
                                            'ORF_seq_similarity',
                                            'NMD_status',
                                            'domains_identified',
                                            #'IDR_identified',
                                            #'IDR_type',
                                            'signal_peptide_identified'),
                                          alpha=0.05,
                                          dIFcutoff=0.1,
                                          onlySigIsoforms=FALSE,
                                          ntCutoff=50,
                                          ntFracCutoff=NULL,
                                          ntJCsimCutoff=0.8,
                                          AaCutoff=10,
                                          AaFracCutoff=0.5,
                                          AaJCsimCutoff=0.9,
                                          removeNonConseqSwitches=TRUE,
                                          showProgress=TRUE
  )

  ### Visual analysis
  # Indiviudal switches

  mostSwitchingGene <- extractTopSwitches( ## options included in plot_mode
    SwitchList,
    n = Inf,#100,
    filterForConsequences=TRUE,
    extractGenes=TRUE,
    alpha=0.05,
    dIFcutoff = 0.1,
    inEachComparison=FALSE, # not included
    sortByQvals=TRUE
  )

  ## This function enables a full analysis of a specific gene containing an isoform switch

  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/gene.pdf', onefile = FALSE, height=6, width = 9)
  switchPlot(
    ### Core arguments
    SwitchList,
    gene = "PRMT5",#mostSwitchingGene$gene_id,
    #isoform_id = NULL,
    condition1 = "cancer", #mostSwitchingGene$condition_1,
    condition2 = "health", #mostSwitchingGene$condition_2,
    ### Advanced arguments
    IFcutoff = 0.05,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    rescaleTranscripts = TRUE,
    reverseMinus = TRUE,
    addErrorbars = TRUE,
    logYaxis = FALSE,
    localTheme = theme_bw(base_size = 8),
  )
  dev.off()


  switchPlotTopSwitches( SwitchList,
                        alpha = 0.05,
                        dIFcutoff = 0.1,
                        onlySigIsoforms = FALSE,
                        n= 10, #inf,
                        sortByQvals=TRUE,
                        filterForConsequences = FALSE,
                        pathToOutput = "/home/laptop/stringtie_test/isoformSwitch_outputs/outputs", #getwd(),
                        splitComparison=FALSE,
                        splitFunctionalConsequences = TRUE,
                        IFcutoff=0.05,
                        fileType = "png",
  )

  pdf(file = '/home/laptop/Galaxy/tools-iuc/tools/isoformswitchanalyzer/test-data/extractConsequencesSummary.pdf', onefile = FALSE, height=6, width = 10)
  consequenceSummary <- extractConsequenceSummary(
    SwitchList,
    consequencesToAnalyze='all',
    includeCombined=FALSE, # not included
    asFractionTotal=FALSE,
    alpha=0.05,
    dIFcutoff=0.1,
    plot=TRUE,
    plotGenes=FALSE,
    simplifyLocation = TRUE,
    removeEmptyConsequences = FALSE,
    returnResult=TRUE, # not included
    localTheme=theme_bw()
  )
  dev.off()


  pdf(file = '/home/laptop/Galaxy/tools-iuc/tools/isoformswitchanalyzer/test-data/consequencesEnrichment.pdf', onefile = FALSE, height=6, width = 9)
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

  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/splicingEnrichment.pdf', onefile = FALSE, height=6, width = 9)
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

  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/splicingSummary.pdf', onefile = FALSE, height=6, width = 9)
  splicingSummary <- extractSplicingSummary(
    SwitchList,
    splicingToAnalyze = 'all',
    asFractionTotal = FALSE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    plot = TRUE,
    plotGenes = FALSE,
    localTheme = theme_bw(),
    returnResult = TRUE
  )
  dev.off()

  ### Volcano like plot:
  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/volcanoPlot.pdf', onefile = FALSE, height=6, width = 9)
  ggplot(data=SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
    geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
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
  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/switchGene.pdf', onefile = FALSE, height=6, width = 9)
  ggplot(data=SwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
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

  pdf(file = '/home/laptop/stringtie_test/isoformSwitch_outputs/outputs/splicingGenomewide.pdf', onefile = FALSE, height=6, width = 9)
  splicingGenomeWide <- extractSplicingGenomeWide(
    SwitchList,
    featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
    splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
    plot=TRUE,
    returnResult=TRUE  # Preventing the summary statistics to be returned as a data.frame
  )
  dev.off()
}