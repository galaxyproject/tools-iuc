# Test 1: Basic FASTA output from BAM
samtools fasta -f 0 -F 2304 -G 0 <(samtools sort -n samtools_fastx-in1.bam) > samtools_fastx-out1.fasta

# Test 2: FASTQ output with r0, r1, r2 splits from BAM
samtools sort -n samtools_fastx-in2.bam | samtools fastq -0 samtools_fastx-out2-1.fastq -1 samtools_fastx-out2-2.fastq -2 samtools_fastx-out2-3.fastq -

# Test 3: FASTA output with r0, r1, r2 splits from SAM
samtools sort -n samtools_fastx-in3.sam | samtools fasta -0 samtools_fastx-out3-1.fasta -1 samtools_fastx-out3-2.fasta -2 samtools_fastx-out3-3.fasta -

# Test 4: Compressed FASTA output from BAM
samtools fasta samtools_fastx-in1.bam | gzip > samtools_fastx-out1.fasta.gz

# Test 5: Compressed FASTQ output with r0, r1, r2 splits from BAM
samtools sort -n samtools_fastx-in2.bam -T tmp | samtools fastq -0 >(gzip > samtools_fastx-out2-1.fastq.gz) -1 >(gzip > samtools_fastx-out2-2.fastq.gz) -2 >(gzip > samtools_fastx-out2-3.fastq.gz) -

# Test 6: Compressed FASTA output with r0, r1, r2 splits from SAM
samtools sort -n samtools_fastx-in3.sam -T tmp | samtools fasta -0 >(gzip > samtools_fastx-out3-1.fasta.gz) -1 >(gzip > samtools_fastx-out3-2.fasta.gz) -2 >(gzip > samtools_fastx-out3-3.fasta.gz) -

# Test 7: Basic 2 output test without singleton tracking
samtools fastq -1 1.1.fq.expected -2 1.2.fq.expected bam2fq.001.sam > 1.stdout.expected

# Test 8: Basic 2 output test with singleton tracking but no singleton
samtools fastq -s 2.s.fq.expected -1 2.1.fq.expected -2 2.2.fq.expected bam2fq.001.sam > 2.stdout.expected

# Test 9: Basic 2 output test with singleton tracking with a singleton in the middle
samtools fastq -s 3.s.fq.expected -1 3.1.fq.expected -2 3.2.fq.expected bam2fq.002.sam > 3.stdout.expected

# Test 10: Basic 2 output test with singleton tracking with a singleton as last read
samtools fastq -s 4.s.fq.expected -1 4.1.fq.expected -2 4.2.fq.expected bam2fq.003.sam > 4.stdout.expected

# Test 11: Tag output test with barcode index
samtools fastq --barcode-tag BC --index-format 'n2i2' --i1 bc.fq.expected -s 4.s.fq.expected -1 4.1.fq.expected -2 4.2.fq.expected bam2fq.004.sam > 4.stdout.expected

# Test 12: Test -O flag with no OQ tags
samtools fastq -O --barcode-tag BC --index-format 'n2i2' --i1 bc.fq.expected -s 4.s.fq.expected -1 4.1.fq.expected -2 4.2.fq.expected bam2fq.004.sam > 4.stdout.expected

# Test 13: Test -O flag with OQ tags
samtools fastq -O --barcode-tag BC --index-format 'n2i2' --i1 bc10.fq.expected -s 10.s.fq.expected -1 10.1.fq.expected -2 10.2.fq.expected bam2fq.010.sam > 2.stdout.expected

# Test 14: Tag output test with separators and -N flag
samtools fastq --barcode-tag BC -N --index-format 'n*i*' --i1 bc_split.fq.expected -s 5.s.fq.expected -1 5.1.fq.expected -2 5.2.fq.expected bam2fq.005.sam > 2.stdout.expected

# Test 15: -t flag
samtools fastq -N -t -s 6.s.fq.expected -1 6.1.fq.expected -2 6.2.fq.expected bam2fq.005.sam > 2.stdout.expected

# Test 16: -T flag
samtools fastq -N -t -T MD,ia -s 7.s.fq.expected -1 7.1.fq.expected -2 7.2.fq.expected bam2fq.005.sam > 2.stdout.expected

# Test 17: -i flag with no index
samtools fastq -N -t -i -T MD,ia --index-format 'n2i2' --i1 11.i.fq.expected --i2 11.i2.fq.expected -s 11.s.fq.expected -1 11.1.fq.expected -2 11.2.fq.expected bam2fq.005.sam > 2.stdout.expected

# Test 18: -i flag with index
samtools fastq --barcode-tag BC -i --index-format 'n2i2' --i1 8.i.fq.expected --i2 8.i2.fq.expected -s 8.s.fq.expected -1 8.1.fq.expected -2 8.2.fq.expected bam2fq.004.sam > 2.stdout.expected