echo "> host" > host.fa
head -c 10000000 /dev/urandom | tr -dc 'ACGT' | head -c 100000 >> host.fa
echo "> contaminant" > contaminant.fa
head -c 10000000 /dev/urandom | tr -dc 'ACGT' | head -c 100000 >> contaminant.fa

art_illumina   -ss HS25   -i host.fa   -p   -l 75   -f 1   -m 200   -s 10   -o host.fq
art_illumina   -ss HS25   -i contaminant.fa   -p   -l 75   -f 1   -m 200   -s 10   -o contaminant.fq

cat host.fq1.fq contaminant.fq1.fq  > host_and_contaminant.fq1.fq
cat host.fq2.fq contaminant.fq2.fq  > host_and_contaminant.fq2.fq

gzip -c host_and_contaminant.fq1.fq > host_and_contaminant.fq1.fq.gz
gzip -c host_and_contaminant.fq2.fq > host_and_contaminant.fq2.fq.gz

# use 10-mers to reduce size
bmtool -d host.fa -o host.bitmask -w 10
srprism mkindex -i host.fa -o host.srprism
makeblastdb -in host.fa -dbtype nucl