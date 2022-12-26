samtools ampliconclip -b eboVir3.1.bed eboVir3.bam | samtools collate -@ 0 -O -u - | samtools fixmate -@ 0 -u - - | samtools sort -o eboVir3.clipped.bam
samtools ampliconclip --strand -b eboVir3.1.bed eboVir3.bam | samtools collate -@ 0 -O -u - | samtools fixmate -@ 0 -u - - | samtools sort -o eboVir3.clipped.strand.bam
samtools ampliconclip --hard-clip -b eboVir3.1.bed eboVir3.bam | samtools collate -@ 0 -O -u - | samtools fixmate -@ 0 -u - - | samtools sort -o eboVir3.hardclipped.bam
samtools ampliconclip --both-ends --no-excluded --strand --fail-len 30 -b eboVir3.1.bed eboVir3.bam | samtools collate -@ 0 -O -u - | samtools fixmate -@ 0 -u - - | samtools sort -o eboVir3.clipped.strand_gt30.bam
