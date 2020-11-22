#!/usr/bin/env R

library(vcfR)
library(tidyverse)

input.vcf <- "variants.vcf"
output.tab <- "variants.tabular"


read.vcf <- read.vcfR(input.vcf)
vcf.fix <- as_tibble(getFIX(read.vcf)) %>%
    select(CHROM, POS, REF, ALT, FILTER)


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
                  'EFF[*].AA', 'EFF[*].GENE')

write.table(vcf.info,
            file=output.tab,
            quote=F, sep='\t', row.names=F)