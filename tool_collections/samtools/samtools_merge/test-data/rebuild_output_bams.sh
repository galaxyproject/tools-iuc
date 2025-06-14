samtools merge -s 1 -f 1.merge.expected.bam <(samtools sort -O sam test_input_1_a.sam) <(samtools sort -O sam test_input_1_b.sam)  <(samtools sort -O sam test_input_1_c.sam)
samtools merge -s 1 -f 2.merge.expected.bam test_input_1_a.bam test_input_1_b.bam test_input_1_c.bam
samtools merge -s 1 -f 3.merge.expected.bam test_input_1_b.bam 
samtools merge -s 1 -f -c -p 4.merge.expected.bam test_input_1_a.bam test_input_1_b.bam
samtools merge -s 1 -f 5.merge.expected.bam test_input_1_a_regex.sam test_input_1_b_regex.sam
samtools merge -s 1 -f -L test_input_1_a.bed 6.merge.expected.bam test_input_1_a.bam test_input_1_b.bam test_input_1_c.bam
