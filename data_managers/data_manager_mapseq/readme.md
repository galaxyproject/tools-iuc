# test data explained

## the tool will expect one file each with the corresponing endings in the DB folder

Trimmed version of: 
* OTU table (id for each taxon) (*.otu)
* Ref. fasta DB (*.fasta)
* Taxan assignemnt of each ref. DB sequence (*.txt)
* clustering of the ref. sequences (starting with 0) corresponding to the ref. sequences (*.mscluster)

## Command to get DBs

```
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_lsu-20200130.tar.gz
mkdir temp
tar xvzf silva_ssu-20200130.tar.gz -C temp
mv temp/* silva_ssu-20200130
```