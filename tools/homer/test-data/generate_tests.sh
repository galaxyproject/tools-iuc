#! /usr/bin/bash
## Generate input data:
if [ ! -e test-data/small_simplified.gtf ]; then
  wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz -O /tmp/Mus_musculus.GRCm38.102.gtf.gz
  zcat /tmp/Mus_musculus.GRCm38.102.gtf.gz | head -n 5 > test-data/small.gtf
  zcat /tmp/Mus_musculus.GRCm38.102.gtf.gz | awk -v OFS="\t" -v start=74667792 -v end=74748393 '$1 == "2" && $5 > start && $4 < end{print "chr"$0}' >> test-data/small.gtf
  # annotatePeaks.pl gives different results all time. I need to simplify the gtf.
  cat test-data/small.gtf | grep -v -P "ENSMUST00000152027|ENSMUST00000156342|ENSMUST00000139005|ENSMUST00000144544|ENSMUST00000111982|ENSMUST00000140666|ENSMUST00000190553|ENSMUST00000132326|ENSMUST00000047830|ENSMUST00000047904|ENSMUST00000111980|ENSMUSG00000065500|ENSMUSG00000100642" > test-data/small_simplified.gtf
fi
if [ ! -e test-data/CTCF_peaks.bed ]; then
  wget https://raw.githubusercontent.com/lldelisle/scriptsForWilleminEtAl2021/main/CTCF/E12_Limbs_Wt_CTCF_colored.bed -O test-data/CTCF_peaks.bed
fi

. <(planemo conda_env homer_gtf_to_annotation.xml)
echo "$(which homer)"
## homer_gtf_to_annotation
## First test
parseGTF.pl test-data/small.gtf ann -features exon start_codon stop_codon > /tmp/annotations.txt
assignGenomeAnnotation /tmp/annotations.txt /tmp/annotations.txt -prioritize test-data/annotations_default.txt > test-data/annotations_default_stats.txt
## Second test
parseGTF.pl test-data/small.gtf ann -features exon start_codon > /tmp/annotations2.txt
assignGenomeAnnotation /tmp/annotations2.txt /tmp/annotations2.txt -prioritize test-data/annotations_exon_start.txt > test-data/annotations_exon_start_stats.txt
## Third test
parseGTF.pl test-data/small.gtf ann -features exon start_codon stop_codon -annTSSstartOffset -50 -annTSSendOffset 50 -annTTSstartOffset -50 -annTTSendOffset 50 > /tmp/annotations3.txt
assignGenomeAnnotation /tmp/annotations3.txt /tmp/annotations3.txt -prioritize test-data/annotations_small_TSSTTS.txt > test-data/annotations_small_TSSTTS_stats.txt

## For annotatePeaks.pl
parseGTF.pl test-data/small_simplified.gtf ann -features exon start_codon stop_codon > /tmp/annotations.txt
assignGenomeAnnotation /tmp/annotations.txt /tmp/annotations.txt -prioritize test-data/annotations_default_simplified.txt

## homer_annotatePeaks
## First test
annotatePeaks.pl test-data/CTCF_peaks.bed none -gtf test-data/small_simplified.gtf -ann test-data/annotations_default_simplified.txt > test-data/CTCF_peaks_first.txt
## Second test
annotatePeaks.pl test-data/CTCF_peaks.bed none -ann test-data/annotations_default.txt > test-data/CTCF_peaks_second.txt
## Third test
annotatePeaks.pl test-data/CTCF_peaks.bed none -gtf test-data/small_simplified.gtf > test-data/CTCF_peaks_third.txt
## Fourth test
annotatePeaks.pl test-data/fake_phix_peaks.bed test-data/phiX174.fasta -CpG > test-data/phiXcpg.txt
## Fifth test
annotatePeaks.pl test-data/fake_phix_peaks.bed none > test-data/phiX_nothing.txt
