#!/usr/bin/env R

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(tidyverse))

tsvEFFfromVCF <- function(input.vcf, output.tab){
    read.vcf <- readVcf(input.vcf)
    chrom.pos <- data.frame(read.vcf@rowRanges)[,c("seqnames", "start")]
    ref.alt.filter <- read.vcf@fixed[,c("REF","ALT","FILTER")]
    dp.af <- read.vcf@info[c("DP","AF")]

    ## Unwrap the DNAStringList
    ref.alt.filter <- data.frame(
        REF = as.character(ref.alt.filter$REF),
        ALT = sapply(1:nrow(ref.alt.filter), function(i){
            as.character(ref.alt.filter$ALT[[i]])
        }),
        FILTER = as.character(ref.alt.filter$FILTER))
    ##
    ## Don't unwrap EFF yet, we need to preserve rows
    eff <- read.vcf@info["EFF"]

    stopifnot(nrow(chrom.pos) == nrow(ref.alt.filter))
    stopifnot(nrow(ref.alt.filter) == nrow(dp.af))
    stopifnot(nrow(dp.af) == nrow(eff))

    ## EFF data contains nested constructs we need to unify all
    ## data sources first, and then explode the EFF column.
    united <- as_tibble(cbind(chrom.pos, ref.alt.filter, dp.af, eff))
    united.explodeRows <- unnest(united, cols = c(EFF))

    united.explodeRows <- united.explodeRows %>%
        dplyr::mutate(CHROM=seqnames, POS=start) %>%
        dplyr::select(CHROM, POS, REF, ALT, FILTER, DP, AF, EFF)

    vcf.info <- united.explodeRows %>%
        separate(EFF, sep="[(|)]",
                 into=c('EFF[*].EFFECT', 'EFF[*].IMPACT', 'EFF[*].FUNCLASS',
                        'codon.change', 'EFF[*].AA', 'AA.length',
                        'EFF[*].GENE', 'trans.biotype', 'gene.coding',
                        'trans.id', 'exon.rank', 'gt.num', 'warnings')) %>%
        dplyr::select('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'DP', 'AF',
                      'EFF[*].EFFECT', 'EFF[*].IMPACT', 'EFF[*].FUNCLASS',
                      'EFF[*].AA', 'EFF[*].GENE') %>%
        ## now we deduplicate any rows that arise from subselecting columns
        dplyr::distinct()

    ## At this point, we would still have rows which share a POS and ALT pair
    ## which could be problematic for the heatmap plot later.
    ##
    ## This is not something to worry about here, and is resolved in the heatmap
    ## script later.

    write.table(vcf.info, file=output.tab, quote=F, sep='\t', row.names=F)
}


                                        # M A I N
stopifnot(exists("samples"))

for (i in 1:nrow(samples)){
    entry = samples[i,];
    if (entry$exts %in% c("vcf", "vcf.gz")){
        in.vcf = entry$files
        out.tsv = tempfile()
        tsvEFFfromVCF(in.vcf, out.tsv)
        ## point to the new file
        samples[i,]$files = out.tsv
        print(paste(entry$ids, ": converted from VCF to tabular"))
    } else {
        print(paste(entry$ids, ": already tabular"))
    }
}
