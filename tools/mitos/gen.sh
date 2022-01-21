## Activate MITOS(1) conda env
conda create -y --quiet --override-channels --channel iuc --channel conda-forge --channel bioconda --channel defaults --name mitos2 mitos=1.0.5 zip

conda activate mitos1

tmp=$(mktemp -d)
runmitos.py --input test-data/NC_012920.fasta --code 2 --outdir $tmp --refdir '/home/maze/workspace/tools-iuc/tools/mitos/test-data/mitos1-refdata/' && zip -9 -y -r $tmp/result.zip $tmp

for i in bed faa fas geneorder gff seq
do
	cp $tmp/result.$i test-data/NC_012920.$i
done
cp $tmp/result test-data/NC_012920.mito
cp $tmp/plots/trnA-*.svg test-data/NC_012920_trnA.svg
cp $tmp/plots/prot.pdf test-data/NC_012920_prot.pdf
cp $tmp/plots/rna.pdf test-data/NC_012920_ncrna.pdf
rm -rf $tmp

conda deactivate

## Activate MITOS2 conda env
conda create -y --quiet --override-channels --channel iuc --channel conda-forge --channel bioconda --channel defaults --name mitos2 mitos=2.0.6 zip

conda activate mitos2

# data for 1st test
tmp=$(mktemp -d)
runmitos.py --input test-data/NC_012920.fasta --code 2 --outdir $tmp --refdir './test-data/' --refseqver 'refseq63m/' --intron 0 --oril 0 --orih 0 --finovl 50  --fragovl 0.2 --fragfac 10.0  --evalue 2.0 --cutoff 0.50 --clipfac 10.0     --ncev 0.01  --maxtrnaovl 50 --maxrrnaovl 50  --noplots
cp $tmp/result.bed test-data/mitos2_NC_012920.bed
rm -rf $tmp

# data for 3rd test
tmp=$(mktemp -d)
runmitos.py --input test-data/NC_012920.fasta --code 2 --outdir $tmp --refdir './test-data/' --refseqver 'refseq63m/' --intron 0 --oril 0 --orih 0 --evalue 3.0 --cutoff 0.49 --clipfac 9.0 --ncbicode --alarab --oldstst  --ncev 0.1 --sensitive --maxtrnaovl 51 --maxrrnaovl 49 && zip -9 -y -r $tmp/result.zip $tmp
for i in faa fas geneorder gff mitos seq
do
	cp $tmp/result.$i test-data/mitos2_NC_012920.$i
done
cp $tmp/plots/prot.pdf test-data/mitos2_NC_012920_prot.pdf
cp $tmp/plots/rna.pdf test-data/mitos2_NC_012920_ncrna.pdf
rm -rf $tmp

conda deactivate
