options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(scPipe)
    library(SingleCellExperiment)
    library(optparse)
    library(readr)
    library(ggplot2)
    library(plotly)
    library(DT)
    library(scater)
    library(scran)
    library(scales)
    library(Rtsne)
})

option_list <- list(
    make_option(c("-fasta","--fasta"), type="character", help="Genome fasta file"),
    make_option(c("-exons","--exons"), type="character", help="Exon annotation gff3 file"),
    make_option(c("-barcodes","--barcodes"), type="character", help="Cell barcodes csv file"),
    make_option(c("-read1","--read1"), type="character", help="Read 1 fastq.gz"),
    make_option(c("-read2","--read2"), type="character", help="Read 2 fastq.gz"),
    make_option(c("-samplename","--samplename"), type="character", help="Name to use for sample"),
    make_option(c("-bs1","--bs1"), type="integer", help="Barcode start in Read 1"),
    make_option(c("-bl1","--bl1"), type="integer", help="Barcode length in Read 1"),
    make_option(c("-bs2","--bs2"), type="integer", help="Barcode start in Read 2"),
    make_option(c("-bl2","--bl2"), type="integer", help="Barcode length in Read 2"),
    make_option(c("-us2","--us2"), type="integer", help="UMI start in Read 2"),
    make_option(c("-ul2","--ul2"), type="integer", help="UMI length in Read 2"),
    make_option(c("-report","--report"), type="logical", help="HTML report of plots"),
    make_option(c("-rdata","--rdata"), type="logical", help="Output RData file")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

fa_fn = args$fasta
anno_fn = args$exons
fq_R1 = args$read1
fq_R2 = args$read2

# Outputs
out_dir = "."
combined_fastq = file.path(out_dir, "combined.fastq")
aligned_bam = file.path(out_dir, "aligned.bam")
mapped_bam = file.path(out_dir, "aligned.mapped.bam")
read_structure = list(
    bs1 = args$bs1,   # barcode start position in fq_R1, -1 indicates no barcode
    bl1 = args$bl1,    # barcode length in fq_R1, 0 since no barcode present
    bs2 = args$bs2,    # barcode start position in fq_R2
    bl2 = args$bl2,   # barcode length in fq_R2
    us = args$us2,    # UMI start position in fq_R2
    ul = args$ul2     # UMI length
)

print("Trimming barcodes")
sc_trim_barcode(combined_fastq,
                fq_R1,
                fq_R2,
                read_structure=read_structure)

print("Building genome index")
Rsubread::buildindex(basename=file.path(out_dir, "fasta_index"), reference=fa_fn)

print("Aligning reads to genome")
Rsubread::align(index=file.path(out_dir, "fasta_index"),
    readfile1=combined_fastq,
    output_file=aligned_bam)

if (!is.null(args$barcodes)) {
  barcode_anno=args$barcodes
} else {
  print("Detecting barcodes")
  # detect 10X barcodes and generate sample_index.csv file
  barcode_anno = "sample_index.csv"
  sc_detect_bc(
      infq=combined_fastq,
      outcsv=barcode_anno, # bacode annotation output file name
      bc_len=read_structure$bl2, # barcode length
      max_reads=5000000,         # only process first 5 million reads
      min_count = 100             # discard cell barcodes with few than 100 hits
  )
}

print("Assigning reads to exons")
sc_exon_mapping(aligned_bam, mapped_bam, anno_fn)

print("De-multiplexing data")
sc_demultiplex(mapped_bam, out_dir, barcode_anno, has_UMI=FALSE)

print("Counting genes")
sc_gene_counting(out_dir, barcode_anno)

print("Creating SingleCellExperiment object")
sce <- create_sce_by_dir(out_dir)

if (!is.null(args$report)) {
print("Creating report")
create_report(sample_name=args$samplename,
   outdir=out_dir,
   r1=fq_R1,
   r2=fq_R2,
   outfq=combined_fastq,
   read_structure=read_structure,
   filter_settings=list(rmlow=TRUE, rmN=TRUE, minq=20, numbq=2),
   align_bam=aligned_bam,
   genome_index=file.path(out_dir, "fasta_index"),
   map_bam=mapped_bam,
   exon_anno=anno_fn,
   stnd=TRUE,
   fix_chr=FALSE,
   barcode_anno=barcode_anno,
   max_mis=1,
   UMI_cor=1,
   gene_fl=FALSE)
}

if (!is.null(args$rdata) ) {
    save(sce, file = file.path(out_dir,"scPipe_analysis.RData"))
}

sessionInfo()
