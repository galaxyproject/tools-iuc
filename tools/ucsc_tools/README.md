# UCSC Genome Browser utilities

Galaxy wrappers for command-line utilities from the UCSC Genome Browser
source tree ([kent](https://github.com/ucscGenomeBrowser/kent)). Each
subdirectory is a self-contained tool with its own `.shed.yml`; there is no
shared `macros.xml`, so each wrapper declares its own `@TOOL_VERSION@`.

Pre-built binaries and the authoritative usage text for every program are at
<https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/FOOTER.txt>.

## Assembly and sequence

| Tool | Tool id | Description |
|---|---|---|
| [`faSplit`](fasplit/) | `fasplit` | Split a FASTA file |
| [`faToTwoBit`](fatotwobit/) | `ucsc_fatotwobit` | convert FASTA to the compact binary 2bit format |
| [`faToVcf`](fatovcf/) | `fatovcf` | Convert a FASTA alignment file to Variant Call Format (VCF) single-nucleotide diffs |
| [`twoBitToFa`](twobittofa/) | `ucsc-twobittofa` | Convert all or part of .2bit file to FASTA |

## Track files for a hub

| Tool | Tool id | Description |
|---|---|---|
| [`bedToBigBed`](bedtobigbed/) | `ucsc_bedtobigbed` | convert a sorted BED to the indexed binary bigBed format |
| [`wigtobigwig`](wigtobigwig/) | `ucsc_wigtobigwig` | bedGraph or Wig to bigWig converter |

## Gene models

| Tool | Tool id | Description |
|---|---|---|
| [`genePredToBed`](genepredtobed/) | `ucsc_genepredtobed` | convert a genePred gene model to BED12 |
| [`gff3ToGenePred`](gff3togenepred/) | `ucsc_gff3togenepred` | convert a GFF3 gene annotation to genePred |

## Hub validation

| Tool | Tool id | Description |
|---|---|---|
| [`hubCheck`](hubcheck/) | `ucsc_hubcheck` | validate a UCSC track hub |

## Chains and nets

| Tool | Tool id | Description |
|---|---|---|
| [`axtChain`](ucsc_axtchain/) | `ucsc_axtchain` | chain together axt or psl alignments |
| [`axtToMaf`](ucsc_axttomaf/) | `ucsc_axtomaf` | Convert dataset from axt to MAF format |
| [`chainAntiRepeat`](ucsc_chainantirepeat/) | `ucsc_chainantirepeat` | Remove repeated chains |
| [`chainNet`](ucsc_chainnet/) | `ucsc_chainnet` | make alignment nets out of alignment chains |
| [`chainPreNet`](ucsc_chainprenet/) | `ucsc_chainprenet` | Remove chains that don't have a chance of being netted |
| [`chainSort`](ucsc_chainsort/) | `ucsc_chainsort` | Sort chains |
| [`chainSwap`](ucsc_chainswap/) | `ucsc_chainswap` | Swap target and query in a chain. |
| [`mafToAxt`](maftoaxt/) | `maftoaxt` | Convert file from MAF to axt format |
| [`netChainSubset`](ucsc_netchainsubset/) | `ucsc_netchainsubset` | Create chain file with subset of chains that appear in the net |
| [`netFilter`](ucsc_netfilter/) | `ucsc_netfilter` | Filter out parts of net |
| [`netSyntenic`](ucsc_netsyntenic/) | `ucsc_netsyntenic` | Add synteny info to a net |
| [`netToAxt`](ucsc_nettoaxt/) | `ucsc_nettoaxt` | Convert net (and chain) to axt format |

## MAF utilities

| Tool | Tool id | Description |
|---|---|---|
| [`mafAddIRows`](maftools/) | `ucsc_mafaddirows` | Add i rows to a MAF file |
| [`mafCoverage`](maftools/) | `ucsc_mafcoverage` | Analyse coverage by MAF files |
| [`mafFetch`](maftools/) | `ucsc_maffetch` | Get overlapping records from an MAF using an index table |
| [`mafFilter`](maftools/) | `ucsc_mafFilter` | Filter MAF files based on various criteria |
| [`mafFrag`](maftools/) | `ucsc_maffrag` | Extract MAF sequences for a region from database |
| [`mafFrags`](maftools/) | `ucsc_maffrags` | Extract MAFs from regions specified in a BED file |
| [`mafGene`](maftools/) | `ucsc_mafgene` | Output protein alignments using maf and genePred |

## Building an assembly hub

A UCSC [track hub](https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html)
serves an assembly and its annotation straight from a web directory, so a
genome the Browser does not host can still be viewed in it. The tools above
cover that path end to end:

1. **Assembly** — `faToTwoBit` packs the genome FASTA into the
   [2bit](https://genome.ucsc.edu/goldenPath/help/twoBit.html) container a
   hub's `genomes.txt` points at with `twoBitPath`. `twoBitToFa` reads it back.
2. **Gene models** — `gff3ToGenePred` converts a GFF3 annotation to
   [genePred](https://genome.ucsc.edu/FAQ/FAQformat.html#format9), and
   `genePredToBed` turns that into BED12.
3. **Tracks** — `bedToBigBed` indexes a **sorted** BED into
   [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html);
   `wigtobigwig` does the same for quantitative data. Both formats are read
   over HTTP byte ranges, so the web server hosting a hub must support them.
4. **Validation** — `hubCheck` fetches a `hub.txt` and reports what the
   Browser would reject. An empty report is a hub that loads.

`bedToBigBed` requires its input to be sorted by chromosome then start; an
unsorted BED is an error, not a warning.

## Chains, nets and MAFs

The remaining tools implement the UCSC pairwise-alignment pipeline, in which
`axtChain` builds chains from alignments, `chainPreNet` / `chainNet` /
`netSyntenic` reduce them to nets, and `netChainSubset` / `netToAxt` /
`axtToMaf` convert the result for downstream multiple alignment.
