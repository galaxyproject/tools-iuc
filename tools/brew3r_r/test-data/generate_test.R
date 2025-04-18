library(GenomicRanges)
input_to_overlap_case1_2_3_4_6_7_8 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
        start = 3,
        end = 25
    ),
    strand = "+",
    gene_id = "geneA",
    transcript_id = "transcriptA",
    type = "exon",
    exon_id = "exonA"
)
big_gr <- NULL
for (i in c(1:5, 7:10)) {
    temp.gr <- input_to_overlap_case1_2_3_4_6_7_8
    temp.gr <- shift(temp.gr, 100 * (i - 1))
    temp.gr$gene_id <- paste0("gene", LETTERS[i])
    temp.gr$transcript_id <- paste0("transcript", LETTERS[i])
    temp.gr$exon_id <- paste0("exon", LETTERS[i])
    temp.gr$exon_number <- 1
    big_gr <- c(big_gr, temp.gr)
}
input_to_overlap_case5_9 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
        start = c(1, 33, 45, 72),
        end = c(25, 40, 60, 75)
    ),
    strand = "+",
    gene_id = "geneA",
    transcript_id = "transcriptA",
    type = "exon",
    exon_id = c("exonA", "exonB", "exonC", "exonD")
)
for (i in c(6, 11)) {
    temp.gr <- input_to_overlap_case5_9
    temp.gr <- shift(temp.gr, 100 * (i - 1))
    temp.gr$gene_id <- paste0("gene", LETTERS[i])
    temp.gr$transcript_id <- paste0("transcript", LETTERS[i])
    temp.gr$exon_id <- paste0("exon", LETTERS[i], letters[1:4])
    temp.gr$exon_number <- 1:4
    big_gr <- c(big_gr, temp.gr)
}
big_gr <- unlist(as(big_gr, "GRangesList"))



input_gr <- c(
    # 1
    GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = c(5, 20),
            end = c(10, 30)
        ),
        strand = "+",
        gene_id = c("gene11", "gene12"),
        transcript_id = c("transcript11", "transcript12"),
        type = "exon",
        exon_id = c("exon11", "exon12")
    ),
    # 2
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 20),
                end = c(10, 25)
            ),
            strand = "+",
            gene_id = c("gene21", "gene22"),
            transcript_id = c("transcript21", "transcript22"),
            type = "exon",
            exon_id = c("exon21", "exon22")
        ),
        100
    ),
    # 3_5
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 20),
                end = c(10, 22)
            ),
            strand = "+",
            gene_id = c("gene31", "gene32"),
            transcript_id = c("transcript31", "transcript32"),
            type = "exon",
            exon_id = c("exon31", "exon32")
        ),
        200
    ),
    # 4
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 5),
                end = c(10, 22)
            ),
            strand = "+",
            gene_id = c("gene41", "gene42"),
            transcript_id = c("transcript41", "transcript42"),
            type = "exon",
            exon_id = c("exon41", "exon42")
        ),
        300
    ),
    # 4bis
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 5),
                end = c(10, 25)
            ),
            strand = "+",
            gene_id = c("gene51", "gene52"),
            transcript_id = c("transcript51", "transcript52"),
            type = "exon",
            exon_id = c("exon51", "exon52")
        ),
        400
    ),
    # 3_5
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 20),
                end = c(10, 22)
            ),
            strand = "+",
            gene_id = c("gene61", "gene62"),
            transcript_id = c("transcript61", "transcript62"),
            type = "exon",
            exon_id = c("exon61", "exon62")
        ),
        500
    ),
    # 6
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(1, 1, 30),
                end = c(10, 10, 40)
            ),
            strand = "+",
            gene_id = "gene71",
            transcript_id = c("transcript71", "transcript72", "transcript72"),
            type = "exon",
            exon_id = c("exon71", "exon71", "exon72")
        ),
        600
    ),
    # 6bis
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(1, 1, 30),
                end = c(10, 10, 40)
            ),
            strand = "+",
            gene_id = c("gene81", "gene82", "gene82"),
            transcript_id = c("transcript81", "transcript82", "transcript82"),
            type = "exon",
            exon_id = c("exon81", "exon82", "exon83")
        ),
        700
    ),
    # 7
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(1, 1, 30),
                end = c(8, 10, 40)
            ),
            strand = "+",
            gene_id = "gene1",
            transcript_id = c("transcript91", "transcript92", "transcript92"),
            type = "exon",
            exon_id = c("exon91", "exon92", "exon93")
        ),
        800
    ),
    # 8
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(1, 1, 30),
                end = c(8, 10, 40)
            ),
            strand = "+",
            gene_id = c("gene101", "gene102", "gene102"),
            transcript_id = c("transcript101", "transcript102", "transcript102"),
            type = "exon",
            exon_id = c("exon101", "exon102", "exon103")
        ),
        900
    ),
    # 9
    shift(
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = c(5, 55),
                end = c(10, 70)
            ),
            strand = "+",
            gene_id = c("gene111", "gene112"),
            transcript_id = c("transcript111", "transcript112"),
            type = "exon",
            exon_id = c("exon111", "exon112")
        ),
        1000
    )
)
## Add convergent genes overlapping a unstranded
input_gr <- c(
    input_gr,
    GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = c(1100, 1110),
            end = c(1105, 1120)
        ),
        strand = c("+", "-"),
        gene_id = c("gene121", "gene122"),
        transcript_id = c("transcript121", "transcript122"),
        type = "exon",
        exon_id = c("exon121", "exon122")
    )
)
big_gr <- c(
    big_gr,
    GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = 1103,
            end = 1113
        ),
        strand = "*",
        gene_id = "geneL",
        transcript_id = "transcriptL",
        type = "exon",
        exon_id = "exonL"
    )
)
input_gr$gene_name <- input_gr$gene_id
input_gr$gene_name[input_gr$gene_id == "gene111"] <- "Gm001"
library(BREW3R.r)
new.gr <- extend_granges(input_gr, big_gr)
library("rtracklayer")
export.gff(input_gr, "input.gtf")
input_gr$exon_id <- NULL
export.gff(input_gr, "input_noexonid.gtf")
export.gff(big_gr, "second_input.gtf")
export.gff(sort(new.gr, ignore.strand = TRUE), "output.gtf")
