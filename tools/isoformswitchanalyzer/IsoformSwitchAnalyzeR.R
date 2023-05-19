# Load the IsoformSwitchAnalyzeR library
library(IsoformSwitchAnalyzeR,
        quietly = TRUE,
        warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)


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
parser <- ArgumentParser(description = "IsoformSwitcheR R script")

parser$add_argument("--modeSelector")
parser$add_argument("--parentDir",  required = FALSE, help = "Parent directory")
parser$add_argument("--readLength",
                    required = FALSE,
                    type = "integer",
                    help = "Read length (required for stringtie)")
parser$add_argument("--annotation", required = FALSE, help = "Annotation")
parser$add_argument("--transcriptome", required = FALSE, help = "Transcriptome")
parser$add_argument(
  "--fixStringTieAnnotationProblem",
  action = "store_true",
  required = FALSE,
  help = "Fix StringTie annotation problem"
)
parser$add_argument("--countFiles", required = FALSE, help = "Count files")
parser$add_argument("--toolSource", required = FALSE, help = "Tool source")
parser$add_argument("--rObject", required = FALSE, help = "R object")
parser$add_argument("--IFcutoff",
                    required = FALSE,
                    type = "numeric",
                    help = "IFcutoff")
parser$add_argument(
  "--geneExpressionCutoff",
  required = FALSE,
  type = "numeric",
  help = "Gene expression cutoff"
)
parser$add_argument(
  "--isoformExpressionCutoff",
  required = FALSE,
  type = "numeric",
  help = "Isoform expression cutoff"
)
parser$add_argument("--alpha",
                    required = FALSE,
                    type = "numeric",
                    help = "")
parser$add_argument("--dIFcutoff",
                    required = FALSE,
                    type = "numeric",
                    help = "dIF cutoff")
parser$add_argument(
  "--onlySigIsoforms",
  required = FALSE,
  action = "store_true",
  help = "Only significative isoforms"
)
parser$add_argument(
  "--filterForConsequences",
  required = FALSE,
  action = "store_true",
  help = "Filter for consequences"
)
parser$add_argument(
  "--removeSingleIsformGenes",
  required = FALSE,
  action = "store_true",
  help = "Remove single isoform genes"
)
parser$add_argument(
  "--keepIsoformInAllConditions",
  required = FALSE,
  action = "store_true",
  help = "Keep isoform in all conditions"
)
parser$add_argument(
  "--correctForConfoundingFactors",
  required = FALSE,
  action = "store_true",
  help = "Correct for confunding factors"
)
parser$add_argument(
  "--overwriteIFvalues",
  required = FALSE,
  action = "store_true",
  help = "Overwrite IF values"
)
parser$add_argument(
  "--reduceToSwitchingGenes",
  required = FALSE,
  action = "store_true",
  help = "Reduce to switching genes"
)
parser$add_argument(
  "--reduceFurtherToGenesWithConsequencePotential",
  required = FALSE,
  action = "store_true",
  help = "Reduce further to genes with consequence potential"
)
parser$add_argument(
  "--keepIsoformInAllConditions2",
  required = FALSE,
  action = "store_true",
  help = "Keep isoform in ll conditions"
)
parser$add_argument("--minORFlength",
                    required = FALSE,
                    type = "integer",
                    help = "")
parser$add_argument("--orfMethod", required = FALSE, help = "ORF methods")
parser$add_argument("--PTCDistance",
                    required = FALSE,
                    type = "integer",
                    help = "")
parser$add_argument(
  "--removeShortAAseq",
  required = FALSE,
  action = "store_true",
  help = "Remove short aminoacid sequences"
)
parser$add_argument(
  "--removeLongAAseq",
  required = FALSE,
  action = "store_true",
  help = "Remove long aminoacid sequences"
)
parser$add_argument(
  "--removeORFwithStop",
  required = FALSE,
  action = "store_true",
  help = "Remove ORF with stop codon"
)
parser$add_argument(
  "--onlySwitchingGenes",
  required = FALSE,
  action = "store_true",
  help = "Only switching genes"
)
parser$add_argument("--analysisMode", required = FALSE, help = "Analyze all isoforms with differential usage or single genes")
parser$add_argument(
  "--genesToPlot",
  type = "integer",
  default = 10,
  required = FALSE,
  help = "Number of genes to plot"
)
parser$add_argument("--gene", required = FALSE, help = "Gene ID to analyze")
parser$add_argument(
  "--sortByQvals",
  action = "store_true",
  required = FALSE,
  help = "Sort genes by Q-val values"
)
parser$add_argument("--countGenes",
                    action = "store_true",
                    required = FALSE,
                    help = "Count genes")
parser$add_argument(
  "--asFractionTotal",
  action = "store_true",
  required = FALSE,
  help = "Plot gene expresson as fraction of total"
)
parser$add_argument("--plotGenes",
                    action = "store_true",
                    required = FALSE,
                    help = "Plot genes instead of isoforms")
parser$add_argument(
  "--simplifyLocation",
  action = "store_true",
  required = FALSE,
  help = "Simplify localtion"
)
parser$add_argument(
  "--removeEmptyConsequences",
  action = "store_true",
  required = FALSE,
  help = "Remove empty consequences"
)
parser$add_argument(
  "--analysisOppositeConsequence",
  action = "store_true",
  required = FALSE,
  help = "Analysi opposite consequences"
)
parser$add_argument("--pathToCPATresultFile",
                    required = FALSE,
                    help = "Path to CPAT result file")
parser$add_argument("--pathToCPC2resultFile",
                    required = FALSE,
                    help = "Path to CPC2 result file")
parser$add_argument("--pathToPFAMresultFile",
                    required = FALSE,
                    help = "Path to PFAM result file")
parser$add_argument("--pathToNetSurfP2resultFile",
                    required = FALSE,
                    help = "Path to NetSurfP2 result file")
parser$add_argument("--pathToSignalPresultFile",
                    required = FALSE,
                    help = "Path to signalP result file")
parser$add_argument("--pathToIUPred2AresultFile",
                    required = FALSE,
                    help = "Path to IUPred2A result file")
parser$add_argument("--codingCutoff",
                    required = FALSE,
                    type = "numeric",
                    help = "Codding cutoff")
parser$add_argument(
  "--removeNoncodingORFs",
  action = "store_true",
  required = FALSE,
  help = "Remove non-coding ORFs"
)
parser$add_argument(
  "--minSignalPeptideProbability",
  required = FALSE,
  type = "numeric",
  help = "Minimul signal peptide probability"
)
parser$add_argument(
  "--smoothingWindowSize",
  type = "integer",
  required = FALSE,
  help = "Smoothing windows size"
)
parser$add_argument(
  "--probabilityCutoff",
  required = FALSE,
  type = "double",
  help = "Probability cutoff"
)
parser$add_argument("--minIdrSize",
                    required = FALSE,
                    type = "integer",
                    help = "Min Idr size")
parser$add_argument(
  "--annotateBindingSites",
  action = "store_true",
  required = FALSE,
  help = "Annotate binding sites"
)
parser$add_argument(
  "--minIdrBindingSize",
  required = FALSE,
  type = "integer",
  help = "Minimun Idr binding size"
)
parser$add_argument(
  "--minIdrBindingOverlapFrac",
  required = FALSE,
  type = "numeric",
  help = ""
)
parser$add_argument("--ntCutoff",
                    required = FALSE,
                    type = "integer",
                    help = "Nucleotide cutoff")
parser$add_argument("--ntFracCutoff",
                    required = FALSE,
                    type = "numeric",
                    help = "Nucleotide fraction cutoff")
parser$add_argument(
  "--ntJCsimCutoff",
  required = FALSE,
  type = "numeric",
  help = "Nucleotide Jaccard simmilarity cutoff"
)
parser$add_argument("--AaCutoff",
                    required = FALSE,
                    type = "integer",
                    help = "Aminoacid cutoff")
parser$add_argument("--AaFracCutoff",
                    required = FALSE,
                    type = "numeric",
                    help = "Aminoacid fraction cutoff")
parser$add_argument(
  "--AaJCsimCutoff",
  required = FALSE,
  type = "numeric",
  help = "Aminoacid Jaccard similarity cutoff"
)
parser$add_argument(
  "--removeNonConseqSwitches",
  action = "store_true",
  required = FALSE,
  help = "Remove switches without consequences"
)
parser$add_argument(
  "--rescaleTranscripts",
  action = "store_true",
  required = FALSE,
  help = "Rescale transcripts"
)
parser$add_argument(
  "--reverseMinus",
  action = "store_true",
  required = FALSE,
  help = "Reverse minus"
)
parser$add_argument(
  "--addErrorbars",
  action = "store_true",
  required = FALSE,
  help = "Add error bars"
)


args <- parser$parse_args()

# Data import
###################

if (args$modeSelector == "data_import") {

  quantificationData <- importIsoformExpression(
    parentDir = args$parentDir,
    addIsofomIdAsColumn = TRUE,
    readLength = args$readLength
  )

  ### Make design matrix
  myDesign <- data.frame(
    sampleID = colnames(quantificationData$abundance)[-1],
    condition = gsub(
      "[[:digit:]]+",
      "",
      colnames(quantificationData$abundance)[-1]
    )
  )

  if (args$toolSource == "stringtie") {
    SwitchList <- importRdata(
      isoformCountMatrix   = quantificationData$counts,
      isoformRepExpression = quantificationData$abundance,
      designMatrix         = myDesign,
      isoformExonAnnoation = args$annotation,
      isoformNtFasta       = args$transcriptome,
      showProgress = TRUE,
      fixStringTieAnnotationProblem = args$fixStringTieAnnotationProblem
    )
  } else {
    SwitchList <- importRdata(
      isoformCountMatrix   = quantificationData$counts,
      isoformRepExpression = quantificationData$abundance,
      designMatrix         = myDesign,
      isoformExonAnnoation = args$annotation,
      isoformNtFasta       = args$transcriptome,
      showProgress = TRUE
    )
  }


  geneCountMatrix <- extractGeneExpression(
    SwitchList,
    extractCounts = TRUE,
    addGeneNames = FALSE,
    addIdsAsColumns = FALSE
  )

  if (args$countFiles == "collection") {

    expressionDF <- data.frame(geneCountMatrix)

    myDesign$condition[length(myDesign$condition)]

    dataframe_factor1 <- expressionDF %>% select(matches(myDesign$condition[1]))
    dataframe_factor2 <- expressionDF %>% select(matches(myDesign$condition[length(myDesign$condition)]))


    lf1 <- as.list(as.data.frame(dataframe_factor1))
    sampleNames1 <- colnames(as.data.frame(dataframe_factor1))

    lf2 <- as.list(as.data.frame(dataframe_factor2))
    sampleNames2 <- colnames(as.data.frame(dataframe_factor2))

    geneNames <- row.names(as.data.frame(expressionDF))


    for (index in seq_along(lf1)) {
      tabular_expression <- data.frame(geneNames, lf1[index])
      colnames(tabular_expression) <-
        c("Geneid", sampleNames1[index])
      filename <-
        paste(sampleNames1[index], "dataset.tabular", sep = "_")
      output_path <- paste("./count_files/factor1/", filename, sep = "")
      write.table(
        tabular_expression,
        output_path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
    for (index in seq_along(lf2)) {
      tabular_expression <- data.frame(geneNames, lf2[index])
      colnames(tabular_expression) <-
        c("Geneid", sampleNames2[index])
      filename <-
        paste(sampleNames2[index], "dataset.tabular", sep = "_")
      output_path <- paste("./count_files/factor2/", filename, sep = "")
      write.table(
        tabular_expression,
        output_path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  } else if (args$countFiles == "matrix") {
    expressionDF <- data.frame(geneCountMatrix)
    geneNames <- row.names(expressionDF)

    expressionDF <- cbind(geneNames, expressionDF)
    write.table(
      as.data.frame(expressionDF),
      "./count_files/matrix.tabular",
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    write.table(
      as.data.frame(myDesign),
      "./count_files/samples.tabular",
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }

  save(SwitchList, file = "SwitchList.Rda")

}

if (args$modeSelector == "first_step") {

  # First part of the analysis
  #############################

  load(file = args$rObject)

  ### Filter
  SwitchList <- preFilter(
    SwitchList,
    IFcutoff = args$IFcutoff,
    geneExpressionCutoff = args$geneExpressionCutoff,
    isoformExpressionCutoff = args$isoformExpressionCutoff,
    removeSingleIsoformGenes = args$removeSingleIsformGenes,
    onlySigIsoforms = args$onlySigIsoforms,
    keepIsoformInAllConditions = args$keepIsoformInAllConditions,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
  )

  ### Test for isoform switches
  SwitchList <- isoformSwitchTestDEXSeq(
    SwitchList,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    correctForConfoundingFactors = args$correctForConfoundingFactors,
    overwriteIFvalues = args$overwriteIFvalues,
    reduceToSwitchingGenes = args$reduceToSwitchingGenes,
    reduceFurtherToGenesWithConsequencePotential = args$reduceFurtherToGenesWithConsequencePotential,
    onlySigIsoforms = args$onlySigIsoforms,
    keepIsoformInAllConditions = args$keepIsoformInAllConditions2,
    showProgress = TRUE,
  )

  SwitchList <- analyzeNovelIsoformORF(
    SwitchList,
    analysisAllIsoformsWithoutORF = TRUE,
    minORFlength = args$minORFlength,
    orfMethod = args$orfMethod,
    PTCDistance = args$PTCDistance,
    startCodons = "ATG",
    stopCodons = c("TAA", "TAG", "TGA"),
    showProgress = TRUE,
  )

  ### Extract Sequences
  SwitchList <- extractSequence(
    SwitchList,
    onlySwitchingGenes = args$onlySwitchingGenes,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    extractNTseq = TRUE,
    extractAAseq = TRUE,
    removeShortAAseq = args$removeShortAAseq,
    removeLongAAseq = args$removeLongAAseq,
    removeORFwithStop = args$removeORFwithStop,
    addToSwitchAnalyzeRlist = TRUE,
    writeToFile = TRUE,
    pathToOutput = getwd(),
    outputPrefix = "isoformSwitchAnalyzeR_isoform",
    forceReExtraction = FALSE,
    quiet = FALSE
  )

  ### Summary
  switchSummary <- extractSwitchSummary(
    SwitchList,
    filterForConsequences = args$filterForConsequences,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    onlySigIsoforms = args$onlySigIsoforms,
  )

  save(SwitchList, file = "SwitchList.Rda")
  write.table(
    switchSummary,
    file = "switchSummary.tsv",
    quote = FALSE,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

}

if (args$modeSelector == "second_step") {

  # Second part of the analysis
  #############################

  load(file = args$rObject)

  ### Add annotation
  if (!is.null(args$pathToCPATresultFile)) {
    SwitchList <- analyzeCPAT(
      SwitchList,
      pathToCPATresultFile = args$pathToCPATresultFile,
      codingCutoff = args$codingCutoff,
      removeNoncodinORFs        = args$removeNoncodingORFs
    )
  }

  if (!is.null(args$pathToCPC2resultFile)) {
    SwitchList <- analyzeCPC2(
      SwitchList,
      pathToCPC2resultFile = args$pathToCPC2resultFile,
      removeNoncodinORFs = args$removeNoncodingORFs
    )
  }

  if (!is.null(args$pathToPFAMresultFile)) {
    pfamFiles <- list.files(path = args$pathToPFAMresultFile,
                            full.names = TRUE)

    SwitchList <- analyzePFAM(SwitchList,
                              pathToPFAMresultFile =  pfamFiles,
                              showProgress = FALSE)
  }

  if (!is.null(args$pathToNetSurfP2resultFile)) {
    netsurfFiles <- list.files(path = args$pathToNetSurfP2resultFile,
                               full.names = TRUE)

    SwitchList <- analyzeNetSurfP2(
      SwitchList,
      pathToNetSurfP2resultFile =  netsurfFiles,
      smoothingWindowSize = args$smoothingWindowSize,
      probabilityCutoff = args$probabilityCutoff,
      minIdrSize = args$minIdrSize,
      showProgress = TRUE
    )
  }

  if (!is.null(args$pathToIUPred2AresultFile)) {
    SwitchList <- analyzeIUPred2A(
      SwitchList,
      pathToIUPred2AresultFile = args$pathToIUPred2AresultFile,
      smoothingWindowSize = args$smoothingWindowSize,
      probabilityCutoff = args$probabilityCutoff,
      minIdrSize = args$minIdrSize,
      annotateBindingSites = args$annotateBindingSites,
      minIdrBindingSize = args$minIdrBindingSize,
      minIdrBindingOverlapFrac = args$minIdrBindingOverlapFrac,
      showProgress = TRUE,
      quiet = FALSE
    )
  }

  if (!is.null(args$pathToSignalPresultFile)) {
    signalpFiles <- list.files(path = args$pathToSignalPresultFile,
                               full.names = TRUE)

    SwitchList <- analyzeSignalP(
      SwitchList,
      pathToSignalPresultFile = signalpFiles,
      minSignalPeptideProbability = args$minSignalPeptideProbability
    )
  }

  SwitchList <- analyzeAlternativeSplicing(
    SwitchList,
    onlySwitchingGenes = args$onlySwitchingGenes,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    showProgress = TRUE
  )

  SwitchList <- analyzeIntronRetention(
    SwitchList,
    onlySwitchingGenes = args$onlySwitchingGenes,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    showProgress = TRUE
  )

  consequences <- c(
    "intron_retention",
    "NMD_status",
    "isoform_seq_similarity",
    "ORF_genomic",
    "tss",
    "tts"
  )

  if (!is.null(args$pathToCPATresultFile) ||
      !is.null(args$pathToCPC2resultFile)) {
    updatedConsequences <- c(consequences, "coding_potential")
    consequences <- updatedConsequences
  }

  if (!is.null(args$pathToPFAMresultFile)) {
    updatedConsequences <- c(consequences, "domains_identified")
    consequences <- updatedConsequences
  }

  if (!is.null(args$pathToSignalPresultFile)) {
    updatedConsequences <- c(consequences, "signal_peptide_identified")
    consequences <- updatedConsequences
  }

  if (!is.null(args$pathToNetSurfP2resultFile) ||
      !is.null(args$pathToIUPred2AresultFile)) {
    updatedConsequences <- c(consequences, "IDR_identified", "IDR_type")
    consequences <- updatedConsequences
  }

  SwitchList <- analyzeSwitchConsequences(
    SwitchList,
    consequencesToAnalyze = consequences,
    alpha = args$alpha,
    dIFcutoff = args$dIFcutoff,
    onlySigIsoforms = args$onlySigIsoforms,
    ntCutoff = args$ntCutoff,
    ntJCsimCutoff = args$ntJCsimCutoff,
    AaCutoff = args$AaCutoff,
    AaFracCutoff = args$AaFracCutoff,
    AaJCsimCutoff = args$AaJCsimCutoff,
    removeNonConseqSwitches = args$removeNonConseqSwitches,
    showProgress = TRUE
  )


  ### Visual analysis
  # Top genes

  if (args$analysisMode == "single") {

    outputFile <- file.path(getwd(), "single_gene.pdf")

    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 9
    )

    switchPlot(
      SwitchList,
      gene = args$gene,
      condition1 = myDesign$condition[1],
      condition2 = myDesign$condition[length(myDesign$condition)],
      IFcutoff = args$IFcutoff,
      dIFcutoff = args$dIFcutoff,
      rescaleTranscripts = args$rescaleTranscripts,
      reverseMinus = args$reverseMinus,
      addErrorbars = args$addErrorbars,
      logYaxis = FALSE,
      localTheme = theme_bw(base_size = 8)
    )
    dev.off()

  } else {
    mostSwitchingGene <-
      extractTopSwitches(
        SwitchList,
        n = Inf,
        filterForConsequences = args$filterForConsequences,
        extractGenes = TRUE,
        alpha = args$alpha,
        dIFcutoff = args$dIFcutoff,
        inEachComparison = FALSE,
        sortByQvals = args$sortByQvals
      )

    write.table(
      mostSwitchingGene,
      file = "mostSwitchingGene.tsv",
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )


    switchPlotTopSwitches(
      SwitchList,
      alpha = args$alpha,
      dIFcutoff = args$dIFcutoff,
      onlySigIsoforms = args$onlySigIsoforms,
      n = args$genesToPlot,
      sortByQvals = args$sortByQvals,
      pathToOutput = getwd(),
      fileType = "pdf"
    )

    outputFile <-
      file.path(getwd(), "extractConsequencesSummary.pdf")

    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 12
    )

    consequenceSummary <- extractConsequenceSummary(
      SwitchList,
      consequencesToAnalyze = "all",
      includeCombined = FALSE,
      asFractionTotal = args$asFractionTotal,
      alpha = args$alpha,
      dIFcutoff = args$dIFcutoff,
      plot = TRUE,
      plotGenes = args$plotGenes,
      simplifyLocation = args$simplifyLocation,
      removeEmptyConsequences = args$removeEmptyConsequences,
      returnResult = TRUE,
      localTheme = theme_bw()
    )
    dev.off()

    write.table(
      consequenceSummary,
      file = "consequencesSummary.tsv",
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )


    outputFile <- file.path(getwd(), "consequencesEnrichment.pdf")
    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 9
    )
    consequenceEnrichment <- extractConsequenceEnrichment(
      SwitchList,
      consequencesToAnalyze = "all",
      alpha = args$alpha,
      dIFcutoff = args$dIFcutoff,
      countGenes = args$countGenes,
      analysisOppositeConsequence = args$analysisOppositeConsequence,
      plot = TRUE,
      localTheme = theme_bw(base_size = 12),
      minEventsForPlotting = 10,
      returnResult = TRUE
    )
    dev.off()

    write.table(
      consequenceEnrichment,
      file = "consequencesEnrichment.tsv",
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )


    outputFile <- file.path(getwd(), "splicingEnrichment.pdf")
    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 9
    )
    splicingEnrichment <- extractSplicingEnrichment(
      SwitchList,
      splicingToAnalyze = "all",
      alpha = args$alpha,
      dIFcutoff = args$dIFcutoff,
      onlySigIsoforms = args$onlySigIsoforms,
      countGenes = args$countGenes,
      plot = TRUE,
      minEventsForPlotting = 10,
      returnResult = TRUE
    )
    dev.off()

    write.table(
      splicingEnrichment,
      file = "splicingEnrichment.tsv",
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )


    outputFile <- file.path(getwd(), "splicingSummary.pdf")
    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 12
    )
    splicingSummary <- extractSplicingSummary(
      SwitchList,
      splicingToAnalyze = "all",
      asFractionTotal = args$asFractionTotal,
      alpha = args$alpha,
      dIFcutoff = args$dIFcutoff,
      onlySigIsoforms = args$onlySigIsoforms,
      plot = TRUE,
      plotGenes = args$plotGenes,
      localTheme = theme_bw(),
      returnResult = TRUE
    )
    dev.off()

    write.table(
      splicingSummary,
      file = "splicingSummary.tsv",
      quote = FALSE,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )


    ### Volcano like plot:
    outputFile <- file.path(getwd(), "volcanoPlot.pdf")

    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 9
    )

    p <- ggplot(data=SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
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
    print(p)
    dev.off()


    ### Switch vs Gene changes:
    outputFile <- file.path(getwd(), "switchGene.pdf")
    pdf(
      file = outputFile,
      height = 6,
      width = 9
    )
    p <- ggplot(data = SwitchList$isoformFeatures,
           aes(x = gene_log2_fold_change, y = dIF)) +
      geom_point(aes(color = abs(dIF) > 0.1 &
                       isoform_switch_q_value < 0.05),
                 size = 1) +
      #facet_wrap(~ condition_2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual("Signficant\nIsoform Switch", values = c("black", "red")) +
      labs(x = "Gene log2 fold change", y = "dIF") +
      theme_bw()
    print(p)
    dev.off()

    outputFile <- file.path(getwd(), "splicingGenomewide.pdf")
    pdf(
      file = outputFile,
      onefile = FALSE,
      height = 6,
      width = 14
    )
    splicingGenomeWide <- extractSplicingGenomeWide(
      SwitchList,
      featureToExtract = "isoformUsage",
      splicingToAnalyze = "all",
      plot = TRUE,
      returnResult = TRUE
    )
    dev.off()
  }
  save(SwitchList, file = "SwitchList.Rda")

}
