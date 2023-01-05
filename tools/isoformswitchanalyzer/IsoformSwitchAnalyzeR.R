# Load the IsoformSwitchAnalyzeR library
library(IsoformSwitchAnalyzeR)


stringtie <- importIsoformExpression(
  parentDir ="/home/laptop/tmporal/",
  addIsofomIdAsColumn = TRUE,
  readLength = 150
)

### Make design matrix
myDesign <- data.frame(
  sampleID = colnames(stringtie$abundance)[-1],
  condition = gsub('[[:digit:]]+', '', colnames(stringtie$abundance)[-1])
)

SwitchList <- importRdata(
  isoformCountMatrix   = stringtie$counts,
  isoformRepExpression = stringtie$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/home/laptop/stringtie_test/annotation_gencode_v42.gtf.gz",
  isoformNtFasta       = "/home/laptop/stringtie_test/gencode.v42.transcripts.fa.gz",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)



geneCountMatrix <- extractGeneExpression(
  SwitchList,
  extractCounts = FALSE, #TRUE,
  addGeneNames = FALSE,
  addIdsAsColumns = FALSE
)

# Prepare for DESeq2 format
l<-as.list(as.data.frame(geneCountMatrix))


# First part of the analysis
#############################

### Filter
SwitchList <- preFilter( SwitchList,
                         IFcutoff=0.01,
                         geneExpressionCutoff = 1,
                         isoformExpressionCutoff = 0,
                         #acceptedGeneBiotype = NULL,
                         #acceptedIsoformClassCode = NULL,
                         removeSingleIsoformGenes = TRUE,
                         #reduceToSwitchingGenes=FALSE,
                         #reduceFurtherToGenesWithConsequencePotential = FALSE,
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

mostSwitchingGene <- extractTopSwitches(
  SwitchList,
  n = Inf,#100,
  filterForConsequences=TRUE,
  extractGenes=TRUE,
  alpha=0.05,
  dIFcutoff = 0.1,
  inEachComparison=FALSE,
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
  includeCombined=FALSE,
  asFractionTotal=FALSE,
  alpha=0.05,
  dIFcutoff=0.1,
  plot=TRUE,
  plotGenes=FALSE,
  simplifyLocation = TRUE,
  removeEmptyConsequences = FALSE,
  returnResult=TRUE,
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
  analysisOppositeConsequence=TRUE,
  plot=TRUE,
  localTheme = theme_bw(base_size = 12),
  minEventsForPlotting = 10,
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
  minEventsForPlotting = 10,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
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






