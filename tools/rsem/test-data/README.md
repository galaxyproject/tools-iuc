# RSEM test fixtures

All sequence and alignment fixtures are synthetic. The genome contains 1,800 bases generated with Python's `random.Random(42)` and successive choices from `ACGT`.

The GTF defines three transcripts on the forward strand. Gene `g1` has transcripts `t1` (600 bases, exons 1–300 and 601–900) and `t2` (450 bases, exons 1–300 and 901–1050). Gene `g2` has transcript `t3` (400 bases, exon 1201–1600). `transcripts.fa` contains the corresponding spliced sequences. `gene_map.tsv` maps the same transcripts to genes.

Each valid alignment fixture represents 120 fragments with 50-base reads. Forty fragments have a unique alignment to `t1`, 30 to `t2`, and 30 to `t3`. Another 20 fragments align equally to the shared first exon of `t1` and `t2`. The expected gene counts are therefore 90 for `g1` and 30 for `g2`.

`paired.bam` contains forward-stranded pairs. `reverse.bam` contains the same fragments with the first mate on the reverse strand. Both files retain adjacent mates and grouped alternative alignments. `single.sam` contains single-end alignments. `genomic.sam` deliberately uses a chromosome reference instead of transcript references. `ungrouped.sam` deliberately separates alignments of the same read.

The expected tables were generated with the Bioconda RSEM 1.3.3 package in `quay.io/biocontainers/rsem:1.3.3--pl5321h077b44d_12`, using one thread. Paired-end runs use the default estimator. The reverse-stranded run also appends names from the GTF. The single-end run uses a fragment-length mean of 150, standard deviation of 20, and an estimated read start position distribution.

RSEM's tagged 1.3.3 source retains the string `RSEM v1.3.1` in `rsem_perl_utils.pm`. Its `--version` output therefore reports that older string, including in the pinned 1.3.3 package.
