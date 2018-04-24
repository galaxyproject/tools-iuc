#!/usr/bin/env bash

# 1st test
mkdir -p 'multiqc_WDir/cutadapt'
cp 'test-data/cutadapt.txt' 'multiqc_WDir/cutadapt/cutadapt.txt'
sed -i.old 's/You are running/This is/' 'multiqc_WDir/cutadapt/cutadapt.txt'

mkdir -p 'multiqc_WDir/fastqc/data_0/file_0'
cp 'test-data/fastqc_1.txt' 'multiqc_WDir/fastqc/data_0/file_0/fastqc_data.txt'
mkdir -p 'multiqc_WDir/fastqc/data_0/file_1'
cp 'test-data/fastqc_2.txt' 'multiqc_WDir/fastqc/data_0/file_1/fastqc_data.txt'

mkdir -p 'multiqc_WDir/flexbar'
cp 'test-data/flexbar.txt' 'multiqc_WDir/flexbar/flexbar.txt'

mkdir -p 'multiqc_WDir/sortmerna'
cp 'test-data/sortmerna.txt' 'multiqc_WDir/sortmerna/sortmerna.txt'

mkdir -p 'multiqc_WDir/trimmomatic'
cp 'test-data/trimmomatic.txt' 'multiqc_WDir/trimmomatic/trimmomatic.txt'

multiqc multiqc_WDir

mv 'multiqc_report.html' 'test-data/pre_alignment_soft_report.html'
mv 'multiqc_data/multiqc.log' 'test-data/pre_alignment_soft_log.txt'
mv 'multiqc_data/multiqc_cutadapt.txt' 'test-data/cutadapt_stats.tabular'
mv 'multiqc_data/multiqc_fastqc.txt' 'test-data/fastqc_stats.tabular'
mv 'multiqc_data/multiqc_flexbar.txt' 'test-data/flexbar_stats.tabular'
mv 'multiqc_data/multiqc_general_stats.txt' 'test-data/pre_alignment_soft_stats.tabular'
mv 'multiqc_data/multiqc_sortmerna.txt' 'test-data/sortmerna_stats.tabular'
mv 'multiqc_data/multiqc_trimmomatic.txt' 'test-data/trimmomatic_stats.tabular'

rm -rf 'multiqc_WDir'
rm -rf 'multiqc_data/'

# 2nd test
mkdir -p 'multiqc_WDir/bismark'
cp 'test-data/bismark.txt' 'multiqc_WDir/bismark/bismark_SE_report.txt'

mkdir -p 'multiqc_WDir/bowtie2'
cp 'test-data/bowtie2_1.txt' 'multiqc_WDir/bowtie2/bowtie2_1.txt'
cp 'test-data/bowtie2_2.txt' 'multiqc_WDir/bowtie2/bowtie2_2.txt'

mkdir -p 'multiqc_WDir/hisat2'
cp 'test-data/hisat2_1.txt' 'multiqc_WDir/hisat2/hisat2_1.txt'
cp 'test-data/hisat2_2.txt' 'multiqc_WDir/hisat2/hisat2_2.txt'

mkdir -p 'multiqc_WDir/hicexplorer'
cp 'test-data/hicexplorer1.log' 'multiqc_WDir/hicexplorer/hicexplorer1'
cp 'test-data/hicexplorer1.log' 'multiqc_WDir/hicexplorer/hicexplorer1_1'
cp 'test-data/hicexplorer2.log' 'multiqc_WDir/hicexplorer/hicexplorer2'

mkdir -p 'multiqc_WDir/kallisto'
cp 'test-data/kallisto_1.txt' 'multiqc_WDir/kallisto/kallisto_1.txt'
cp 'test-data/kallisto_2.txt' 'multiqc_WDir/kallisto/kallisto_2.txt'

#mkdir -p 'multiqc_WDir/salmon/fld_0/file_0'
#cp 'test-data/salmon_flenDist.txt' 'multiqc_WDir/salmon/fld_0/file_0/flenDist.txt'
#mkdir -p 'multiqc_WDir/salmon/fld_1/file_0'
#cp 'test-data/salmon_meta_info.json' 'multiqc_WDir/salmon/fld_1/file_0/meta_info.json'

mkdir -p 'multiqc_WDir/star/log_0'
cp 'test-data/star_log.txt' 'multiqc_WDir/star/log_0/star_log_Log.final.out'
mkdir -p 'multiqc_WDir/star/genecounts_1'
cp 'test-data/star_counts.txt' 'multiqc_WDir/star/genecounts_1/star_counts_ReadsPerGene.out.tab'

mkdir -p 'multiqc_WDir/tophat'
cp 'test-data/tophat.txt' 'multiqc_WDir/tophat/tophat_align_summary.txt'

multiqc multiqc_WDir

mv 'multiqc_report.html' 'test-data/aligner_soft_report.html'
mv 'multiqc_data/multiqc_bismark_alignment.txt' 'test-data/bismark_stats.tabular'
mv 'multiqc_data/multiqc_bowtie2.txt' 'test-data/bowtie2_stats.tabular'
mv 'multiqc_data/multiqc_general_stats.txt' 'test-data/aligner_soft_stats.tabular'
mv 'multiqc_data/multiqc_hisat2.txt' 'test-data/hisat2_stats.tabular'
mv 'multiqc_data/multiqc_kallisto.txt' 'test-data/kallisto_stats.tabular'
mv 'multiqc_data/multiqc_star.txt' 'test-data/star_stats.tabular'
mv 'multiqc_data/multiqc_tophat.txt.txt' 'test-data/tophat_stats.tabular'

rm -rf 'multiqc_WDir'
rm -rf 'multiqc_data/'

# 3rd test
mkdir -p 'multiqc_WDir/bamtools'
cp 'test-data/bamtools.txt' 'multiqc_WDir/bamtools/bamtools.txt'

mkdir -p 'multiqc_WDir/bcftools'
cp 'test-data/bcftools.txt' 'multiqc_WDir/bcftools/bcftools.txt'

mkdir -p 'multiqc_WDir/busco'
cp 'test-data/busco.txt' 'multiqc_WDir/busco/short_summary_busco.txt'

mkdir -p 'multiqc_WDir/deeptools'
cp 'test-data/deeptools_bamPEFragmentSize.txt' 'multiqc_WDir/deeptools/'
cp 'test-data/deeptools_estimateReadFiltering.txt' 'multiqc_WDir/deeptools/'
cp 'test-data/deeptools_plotCoverageOutRawCounts.txt' 'multiqc_WDir/deeptools/'
cp 'test-data/deeptools_plotCoverageStdout.txt' 'multiqc_WDir/deeptools/'
cp 'test-data/deeptools_plotEnrichment.txt' 'multiqc_WDir/deeptools/'
cp 'test-data/deeptools_plotFingerprintOutRawCounts.txt' 'multiqc_WDir/deeptools/'

mkdir -p 'multiqc_WDir/featureCounts'
cp 'test-data/featureCounts.txt' 'multiqc_WDir/featureCounts/featureCounts.summary'

mkdir -p 'multiqc_WDir/gatk'
cp 'test-data/gatk_BaseRecalibrator.txt' 'multiqc_WDir/gatk/gatk_BaseRecalibrator.txt'
cp 'test-data/gatk_varianteval.txt' 'multiqc_WDir/gatk/gatk_varianteval.txt'

mkdir -p 'multiqc_WDir/htseq'
cp 'test-data/htseq.txt' 'multiqc_WDir/htseq/htseq.txt'

mkdir -p 'multiqc_WDir/picard'
cp 'test-data/picard_collectGcBias.txt' 'multiqc_WDir/picard/picard_collectGcBias.txt'
cp 'test-data/picard_CollectInsertSizeMetrics.txt' 'multiqc_WDir/picard/picard_CollectInsertSizeMetrics.txt'
cp 'test-data/picard_MarkDuplicates.txt' 'multiqc_WDir/picard/picard_MarkDuplicates.txt'
cp 'test-data/picard_CollectBaseDistributionByCycle.txt' 'multiqc_WDir/picard/picard_CollectBaseDistributionByCycle.txt'
cp 'test-data/picard_CollectRnaSeqMetrics.txt' 'multiqc_WDir/picard/picard_CollectRnaSeqMetrics.txt'
cp 'test-data/picard_CollectAlignmentSummaryMetrics.txt' 'multiqc_WDir/picard/picard_CollectAlignmentSummaryMetrics.txt'

mkdir -p 'multiqc_WDir/prokka'
cp 'test-data/prokka_1.txt' 'multiqc_WDir/prokka/prokka_1.txt'
cp 'test-data/prokka_2.txt' 'multiqc_WDir/prokka/prokka_2.txt'

mkdir -p 'multiqc_WDir/quast/file_0'
cp 'test-data/quast.tsv' 'multiqc_WDir/quast/file_0/report.tsv'

#mkdir -p 'multiqc_WDir/rsem'
#cp 'test-data/rsem.txt' 'multiqc_WDir/rsem/rsem.cnt'

mkdir -p 'multiqc_WDir/rseqc/read_gc_0'
cp 'test-data/rseqc.txt' 'multiqc_WDir/rseqc/read_gc_0/rseq.GC.xls'

mkdir -p 'multiqc_WDir/samblaster'
cp 'test-data/samblaster.txt' 'multiqc_WDir/samblaster/samblaster.txt'

mkdir -p 'multiqc_WDir/samtools/stats_0'
cp 'test-data/samtools_stats.txt' 'multiqc_WDir/samtools/stats_0/samtools_stats.txt'
mkdir -p 'multiqc_WDir/samtools/flagstat_1'
cp 'test-data/samtools_flagstat.txt' 'multiqc_WDir/samtools/flagstat_1/samtools_flagstat.txt'
mkdir -p 'multiqc_WDir/samtools/idxstats_2'
cp 'test-data/samtools_flagstat.txt' 'multiqc_WDir/samtools/idxstats_2/samtools_idxstats_idxstat'

mkdir -p 'multiqc_WDir/snpeff'
cp 'test-data/snpeff.csv' 'multiqc_WDir/snpeff/snpeff.txt'

mkdir -p 'multiqc_WDir/vcftools/tstv_by_qual_0'
cp 'test-data/vcftools.txt' 'multiqc_WDir/vcftools/tstv_by_qual_0/vcftools.TsTv.qual'

multiqc multiqc_WDir

mv 'multiqc_report.html' 'test-data/post_aligner_soft_report.html'
mv 'multiqc_data/multiqc_bamtools_stats.txt' 'test-data/bamtools_stats.tabular'
mv 'multiqc_data/multiqc_bcftools_stats.txt' 'test-data/bcftools_stats.tabular'
mv 'multiqc_data/multiqc_busco.txt' 'test-data/busco_stats.tabular'
mv 'multiqc_data/multiqc_featureCounts.txt' 'test-data/featureCounts_stats.tabular'
mv 'multiqc_data/multiqc_gatk_varianteval.txt' 'test-data/gatk_varianteval_stats.tabular'
mv 'multiqc_data/multiqc_general_stats.txt' 'test-data/post_aligner_soft_stats.tabular'
mv 'multiqc_data/multiqc_htseq.txt' 'test-data/htseq_stats.tabular'
mv 'multiqc_data/multiqc_picard_AlignmentSummaryMetrics.txt' 'test-data/picard_AlignmentSummaryMetrics_stats.tabular'
mv 'multiqc_data/multiqc_picard_RnaSeqMetrics.txt' 'test-data/picard_RnaSeqMetrics_stats.tabular'
mv 'multiqc_data/multiqc_picard_baseContent.txt' 'test-data/picard_baseContent_stats.tabular'
mv 'multiqc_data/multiqc_picard_dups.txt' 'test-data/picard_dups_stats.tabular'
mv 'multiqc_data/multiqc_picard_insertSize.txt' 'test-data/picard_insertSize_stats.tabular'
mv 'multiqc_data/multiqc_prokka.txt' 'test-data/prokka_stats.tabular'
mv 'multiqc_data/multiqc_quast.txt' 'test-data/quast_stats.tabular'
#mv 'multiqc_data/multiqc_rsem.txt' 'test-data/rsem_stats.tabular'
mv 'multiqc_data/multiqc_samblaster.txt' 'test-data/samblaster_stats.tabular'
mv 'multiqc_data/multiqc_samtools_flagstat.txt' 'test-data/samtools_flagstat_stats.tabular'
mv 'multiqc_data/multiqc_samtools_stats.txt' 'test-data/samtools_stats_stats.tabular'
mv 'multiqc_data/multiqc_snpeff.txt' 'test-data/snpeff_stats.tabular'

rm -rf 'multiqc_WDir'
rm -rf 'multiqc_data/'

# 4th test
mkdir -p 'multiqc_WDir/custom_content_0'
cp 'test-data/cc_ko15.bpc.tab' 'multiqc_WDir/custom_content_0/file_0_0.tsv'
cp 'test-data/cc_wt15.bpc.tab' 'multiqc_WDir/custom_content_0/file_0_1.tsv'

echo "custom_data:" > 'config_file'
echo "    section_0:" >> 'config_file'
echo "        file_format: 'tsv'" >> 'config_file'
echo "        section_name: 'BPC'" >> 'config_file'
echo "        title: 'Base peak chromatogram'" >> 'config_file'
echo "        description: 'Sum of intensity (Y) of the most intense peaks at each retention time(X)'" >> 'config_file'
echo "        plot_type: 'linegraph'" >> 'config_file'
echo "        pconfig:" >> 'config_file'
echo "            id: 'section_0_linegraph'" >> 'config_file'
echo "            ylab: 'Base Peak Intensity'" >> 'config_file'
echo "            xlab: 'Retention Time'" >> 'config_file'
echo "sp:" >> 'config_file'
echo "    section_0:" >> 'config_file'
echo "        fn: 'file_0_*'" >> 'config_file'

multiqc multiqc_WDir -c 'config_file'

mv 'multiqc_report.html' 'test-data/report_manual_custom_content.html'

rm 'config_file'
rm -rf 'multiqc_WDir'
rm -rf 'multiqc_data/'
