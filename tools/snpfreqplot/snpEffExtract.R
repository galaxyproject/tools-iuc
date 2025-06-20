#!/usr/bin/env R

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(tidyverse))

tsv_eff_from_vcf <- function(input_vcf, output_tab) {
    read_vcf <- readVcf(input_vcf) # nolint
    if (!nrow(read_vcf@fixed)) {
        # no variants in file -> just write a valid header line
        write(c("CHROM", "POS", "REF", "ALT", "AF", "EFF[*].GENE", "EFF[*].EFFECT"),
            ncolumns = 7, file = output_tab, sep = "\t"
        )
        return()
    }
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
        FILTER = as.character(ref_alt_filter$FILTER)
    )
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

    ## EFF columns are defined here:
    ## https://pcingola.github.io/SnpEff/se_inputoutput/
    options(warn = -1) ## suppress warnings
    seperated_info <- united_exploderows %>%
        separate(EFF,
            sep = "[(|)]",
            extra = "merge", ## extra values merged into "extra" column
            into = c(
                "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS",
                "codon.change", "EFF[*].AA", "AA.length",
                "EFF[*].GENE", "trans.biotype", "gene.coding",
                "trans.id", "exon.rank", "gt.num", "warnings",
                "extra"
            )
        )
    options(warn = 0)
    ## If there is data that has been dropped or filled-in, we will see it in
    ## the "extra" column if it isn't NA or an empty quote.
    test_missing <- seperated_info %>%
        dplyr::select("CHROM", "POS", "extra") %>%
        replace_na(list(extra = "")) %>%
        filter(extra != "")

    if (nrow(test_missing) > 0) {
        print(test_missing)
        stop("Extra values were not parsed")
    }

    vcf_info <- seperated_info %>%
        dplyr::select(
            "CHROM", "POS", "REF", "ALT", "FILTER", "DP", "AF",
            "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS",
            "EFF[*].AA", "EFF[*].GENE"
        ) %>%
        ## now we de-duplicate any rows that arise from subselecting columns
        dplyr::distinct()

    ## At this point, we would still have rows which share a POS and ALT pair
    ## which could be problematic for the heatmap plot later.
    ##
    ## This is not something to worry about here, and is resolved in the heatmap
    ## script later.
    write.table(vcf_info,
        file = output_tab,
        quote = F, sep = "\t", row.names = F
    )
}


# M A I N
stopifnot(exists("samples"))

for (i in seq_len(nrow(samples))) {
    entry <- samples[i, ]
    if (entry$exts %in% c("vcf", "vcf.gz")) {
        in_vcf <- entry$files
        out_tsv <- paste0(entry$ids, ".tsv") ## use local dir
        tsv_eff_from_vcf(in_vcf, out_tsv)
        ## point to the new file
        samples[i, ]$files <- out_tsv
        message(paste(entry$ids, ": converted from VCF to tabular"))
    } else {
        message(paste(entry$ids, ": already tabular"))
    }
}
