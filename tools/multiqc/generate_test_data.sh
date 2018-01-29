#!/usr/bin/env bash

# 1st test
mkdir multiqc_WDir

mkdir 'multiqc_WDir/cutadapt_0'
cp 'test-data/cutadapt.txt' 'multiqc_WDir/cutadapt_0/cutadapt.txt'
sed -i.old 's/You are running/This is/' 'multiqc_WDir/cutadapt_0/cutadapt.txt'

mkdir -p 'multiqc_WDir/fastqc_1/data_0/file_0'
cp 'test-data/fastqc_1.txt' 'multiqc_WDir/fastqc_1/data_0/file_0/fastqc_data.txt'
mkdir 'multiqc_WDir/fastqc_1/data_0/file_1'
cp 'test-data/fastqc_2.txt' 'multiqc_WDir/fastqc_1/data_0/file_1/fastqc_data.txt'

mkdir 'multiqc_WDir/flexbar_2'
cp 'test-data/flexbar.txt' 'multiqc_WDir/flexbar_2/flexbar.txt'

mkdir 'multiqc_WDir/sortmerna_3'
cp 'test-data/sortmerna.txt' 'multiqc_WDir/sortmerna_3/sortmerna.txt'

mkdir 'multiqc_WDir/trimmomatic_4'
cp 'test-data/trimmomatic.txt' 'multiqc_WDir/trimmomatic_4/trimmomatic.txt'

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
mkdir multiqc_WDir

mkdir 'multiqc_WDir/bismark_0'
cp 'test-data/bismark.txt' 'multiqc_WDir/bismark_0/bismark_SE_report.txt'

mkdir 'multiqc_WDir/bowtie2_1'
cp 'test-data/bowtie2_1.txt' 'multiqc_WDir/bowtie2_1/bowtie2_1.txt'
cp 'test-data/bowtie2_2.txt' 'multiqc_WDir/bowtie2_1/bowtie2_2.txt'

mkdir 'multiqc_WDir/hisat2_3'
cp 'test-data/hisat2_1.txt' 'multiqc_WDir/hisat2_3/hisat2_1.txt'
cp 'test-data/hisat2_2.txt' 'multiqc_WDir/hisat2_3/hisat2_2.txt'

mkdir 'multiqc_WDir/kallisto_4'
cp 'test-data/kallisto_1.txt' 'multiqc_WDir/kallisto_4/kallisto_1.txt'
cp 'test-data/kallisto_2.txt' 'multiqc_WDir/kallisto_4/kallisto_2.txt'

mkdir -p 'multiqc_WDir/salmon_5/fld_0/file_0'
cp 'test-data/salmon.txt' 'multiqc_WDir/salmon_5/fld_0/file_0/flenDist.txt'

mkdir -p 'multiqc_WDir/star_6/log_0'
cp 'test-data/star_log.txt' 'multiqc_WDir/star_6/log_0/star_log_Log.final.out'
mkdir 'multiqc_WDir/star_6/genecounts_1'
cp 'test-data/star_counts.txt' 'multiqc_WDir/star_6/genecounts_1/star_counts_ReadsPerGene.out.tab'

mkdir 'multiqc_WDir/tophat_7'
cp 'test-data/tophat.txt' 'multiqc_WDir/tophat_7/tophat_align_summary.txt'

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
mkdir multiqc_WDir

mkdir 'multiqc_WDir/bamtools_0'
cp 'test-data/bamtools.txt' 'multiqc_WDir/bamtools_0/bamtools.txt'

mkdir 'multiqc_WDir/bcftools_1'
cp 'test-data/bcftools.txt' 'multiqc_WDir/bcftools_1/bcftools.txt'

mkdir 'multiqc_WDir/busco_2'
cp 'test-data/busco.txt' 'multiqc_WDir/busco_2/short_summary_busco.txt'

mkdir 'multiqc_WDir/featureCounts_3'
cp 'test-data/featureCounts.txt' 'multiqc_WDir/featureCounts_3/featureCounts.summary'

mkdir 'multiqc_WDir/gatk_4'
cp 'test-data/gatk_BaseRecalibrator.txt' 'multiqc_WDir/gatk_4/gatk_BaseRecalibrator.txt'
cp 'test-data/gatk_varianteval.txt' 'multiqc_WDir/gatk_4/gatk_varianteval.txt'

mkdir 'multiqc_WDir/htseq_5'
cp 'test-data/htseq.txt' 'multiqc_WDir/htseq_5/htseq.txt'

mkdir 'multiqc_WDir/picard_6'
cp 'test-data/picard_collectGcBias.txt' 'multiqc_WDir/picard_6/picard_collectGcBias.txt'
cp 'test-data/picard_CollectInsertSizeMetrics.txt' 'multiqc_WDir/picard_6/picard_CollectInsertSizeMetrics.txt'
cp 'test-data/picard_MarkDuplicates.txt' 'multiqc_WDir/picard_6/picard_MarkDuplicates.txt'
cp 'test-data/picard_CollectBaseDistributionByCycle.txt' 'multiqc_WDir/picard_6/picard_CollectBaseDistributionByCycle.txt'
cp 'test-data/picard_CollectRnaSeqMetrics.txt' 'multiqc_WDir/picard_6/picard_CollectRnaSeqMetrics.txt'
cp 'test-data/picard_CollectAlignmentSummaryMetrics.txt' 'multiqc_WDir/picard_6/picard_CollectAlignmentSummaryMetrics.txt'

mkdir 'multiqc_WDir/prokka_7'
cp 'test-data/prokka_1.txt' 'multiqc_WDir/prokka_7/prokka_1.txt'
cp 'test-data/prokka_2.txt' 'multiqc_WDir/prokka_7/prokka_2.txt'

mkdir -p 'multiqc_WDir/quast_8/file_0'
cp 'test-data/quast.tsv' 'multiqc_WDir/quast_8/file_0/report.tsv'

#mkdir 'multiqc_WDir/rsem_9'
#cp 'test-data/rsem.txt' 'multiqc_WDir/rsem_9/rsem.cnt'

mkdir -p 'multiqc_WDir/rseqc_10/read_gc_0'
cp 'test-data/rseqc.txt' 'multiqc_WDir/rseqc_10/read_gc_0/rseq.GC.xls'

mkdir 'multiqc_WDir/samblaster_11'
cp 'test-data/samblaster.txt' 'multiqc_WDir/samblaster_11/samblaster.txt'

mkdir -p 'multiqc_WDir/samtools_12/stats_0'
cp 'test-data/samtools_stats.txt' 'multiqc_WDir/samtools_12/stats_0/samtools_stats.txt'
mkdir 'multiqc_WDir/samtools_12/flagstat_1'
cp 'test-data/samtools_flagstat.txt' 'multiqc_WDir/samtools_12/flagstat_1/samtools_flagstat.txt'
mkdir 'multiqc_WDir/samtools_12/idxstats_2'
cp 'test-data/samtools_flagstat.txt' 'multiqc_WDir/samtools_12/idxstats_2/samtools_idxstats_idxstat'

#mkdir 'multiqc_WDir/snpeff_13'
#cp 'test-data/snpeff.csv' 'multiqc_WDir/snpeff_13/snpeff.txt'

mkdir -p 'multiqc_WDir/vcftools_14/tstv_by_qual_0'
cp 'test-data/vcftools.txt' 'multiqc_WDir/vcftools_14/tstv_by_qual_0/vcftools.TsTv.qual'

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
mv 'multiqc_data/multiqc_rseqc.txt' 'test-data/rseqc_stats.tabular'
mv 'multiqc_data/multiqc_samblaster.txt' 'test-data/samblaster_stats.tabular'
mv 'multiqc_data/multiqc_samtools_flagstat.txt' 'test-data/samtools_flagstat_stats.tabular'
mv 'multiqc_data/multiqc_samtools_stats.txt' 'test-data/samtools_stats_stats.tabular'
#mv 'multiqc_data/multiqc_snpeff.txt' 'test-data/snpeff_stats.tabular'

rm -rf 'multiqc_WDir'
rm -rf 'multiqc_data/'

# 4th test
mkdir multiqc_WDir

mkdir 'multiqc_WDir/custom_content_0'
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
