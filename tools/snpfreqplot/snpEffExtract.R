#!/usr/bin/env R

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(tidyverse))

tsv_eff_from_vcf <- function(input_vcf, output_tab) {
    read_vcf <- readVcf(input_vcf)  # nolint
    chrom_pos <- data.frame(read_vcf@rowRanges)[, c("seqnames", "start")]
    ref_alt_filter <- read_vcf@fixed[, c("REF", "ALT", "FILTER")]
    dp_af <- read_vcf@info[c("DP", "AF")]

    ## Unwrap the DNAStringList
    # nolint start
    ref_alt_filter <- data.frame(
        REF = as.character(ref_alt_filter$REF),
        ALT = sapply(seq_len(nrow(ref_alt_filter)), function(i) {
            as.character(ref_alt_filter$ALT[[i]])
        }),
        FILTER = as.character(ref_alt_filter$FILTER))
    # nolint end
    ##
    ## Don't unwrap EFF yet, we need to preserve rows
    eff <- read_vcf@info["EFF"]

    stopifnot(nrow(chrom_pos) == nrow(ref_alt_filter))
    stopifnot(nrow(ref_alt_filter) == nrow(dp_af))
    stopifnot(nrow(dp_af) == nrow(eff))

    ## EFF data contains nested constructs we need to unify all
    ## data sources first, and then explode the EFF column.
    united <- as_tibble(cbind(chrom_pos, ref_alt_filter, dp_af, eff)) # nolint
    united_exploderows <- unnest(united, cols = c(EFF)) # nolint

    united_exploderows <- united_exploderows %>%
        dplyr::mutate(CHROM = seqnames, POS = start) %>%
        dplyr::select(CHROM, POS, REF, ALT, FILTER, DP, AF, EFF)

    vcf_info <- united_exploderows %>%
        separate(EFF, sep = "[(|)]",
                 into = c("EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS",
                        "codon.change", "EFF[*].AA", "AA.length",
                        "EFF[*].GENE", "trans.biotype", "gene.coding",
                        "trans.id", "exon.rank", "gt.num", "warnings")) %>%
        dplyr::select("CHROM", "POS", "REF", "ALT", "FILTER", "DP", "AF",
                      "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS",
                      "EFF[*].AA", "EFF[*].GENE") %>%
        ## now we deduplicate any rows that arise from subselecting columns
        dplyr::distinct()

    ## At this point, we would still have rows which share a POS and ALT pair
    ## which could be problematic for the heatmap plot later.
    ##
    ## This is not something to worry about here, and is resolved in the heatmap
    ## script later.

    write.table(vcf_info, file = output_tab,
                quote = F, sep = "\t", row.names = F)
}


                                        # M A I N
stopifnot(exists("samples"))

for (i in seq_len(nrow(samples))) {
    entry <- samples[i, ];
    if (entry$exts %in% c("vcf", "vcf.gz")) {
        in_vcf <- entry$files
        out_tsv <- tempfile(pattern=entry$ids, fileext=".tsv")
        tsv_eff_from_vcf(in_vcf, out_tsv)
        ## point to the new file
        samples[i, ]$files <- out_tsv
        print(paste(entry$ids, ": converted from VCF to tabular"))
    } else {
        print(paste(entry$ids, ": already tabular"))
    }
}
