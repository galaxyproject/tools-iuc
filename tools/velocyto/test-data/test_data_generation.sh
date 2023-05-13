cd /home/ldelisle/Documents/mygit/tools-iuc/tools/rgrnastar/test-data
STAR --runMode genomeGenerate --genomeDir 'tempstargenomedir' --genomeFastaFiles filtered3.Homo_sapiens.GRCh38.dna.chromosome.21.fa --sjdbOverhang 100 --sjdbGTFfile filtered3.Homo_sapiens.GRCh38.100.chr21.gtf --genomeSAindexNbases 4
STAR --genomeLoad NoSharedMemory --genomeDir tempstargenomedir --soloType CB_UMI_Simple --readFilesIn pbmc_1k_v2_L001.R2.10k.fastq.gz pbmc_1k_v2_L001.R1.10k.fastq.gz --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloCBwhitelist filtered.barcodes.txt --soloBarcodeReadLength 1 --readFilesCommand zcat --outSAMattributes NH AS HI nM CB UB --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ../../velocyto/test-data/STARsolo_allSAMat.bam
mv Solo.out/Gene/filtered/barcodes.tsv ../../velocyto/test-data/

cd ../../velocyto/test-data/
mkdir -p sample/outs/filtered_gene_bc_matrices/whatever/
ln -s $PWD/barcodes.tsv sample/outs/filtered_gene_bc_matrices/whatever/
ln -s $PWD/STARsolo_allSAMat.bam sample/outs/possorted_genome_bam.bam
velocyto run10x sample ./filtered3.Homo_sapiens.GRCh38.100.chr21.gtf
rm -r sample/velocyto
velocyto run -c -U ./STARsolo_allSAMat.bam ./filtered3.Homo_sapiens.GRCh38.100.chr21.gtf
velocyto run -c -U ./STARsolo_allSAMat.bam ./STARsolo_allSAMat.bam ./filtered3.Homo_sapiens.GRCh38.100.chr21.gtf
