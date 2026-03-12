#/bin/bash

set -e

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz


# create blast DB
# diamond expects 1234 in the tax data: https://github.com/bbuchfink/diamond/blob/56214dfcb4278f08e935147e8dbea7672997386e/src/data/blastdb/blastdb.cpp#L170
# more precisely in the taxdb.bt* files (which are here constructed from the dmp files)
# we also add the path to the root (guess not needed strictly)
# ideally 1234 should also in the sqlite DB, but since the taxon is not in the fasta it should be fine
sqlite3 taxonomy4blast.sqlite3 "SELECT * FROM TaxidInfo;" | sed 's/|/\n/g' | sort -n -u | sed 's/^/^/; s/$/\\s/' > grep.txt
echo "^1234\\s" >> grep.txt
echo "^189779\\s" >> grep.txt
echo "^189778\\s" >> grep.txt
echo "^203693\\s" >> grep.txt
echo "^40117\\s" >> grep.txt
echo "^3379134\\s" >> grep.txt
echo "^2\\s" >> grep.txt

grep -f grep.txt names.dmp > ../ncbi_taxonomy/names.dmp
grep -f grep.txt nodes.dmp > ../ncbi_taxonomy/nodes.dmp

python taxdb.py 
makeblastdb -in db.fasta -parse_seqids -blastdb_version 5 -taxid_map map.txt -title "cox1 blastp DB" -dbtype prot

# create small dmnd data base with taxonomy
# the important thing to get a small DB is to have consecutive taxIDs
# NOTE: filter_and_map_ids modifies taxIDs (to get a small file), i.e. taxIDs will be different from tests using BLAST DB from above
python filter_and_map_ids.py names.dmp nodes.dmp prot.accession2taxid ../names.dmp ../nodes.dmp ../prot.accession2taxid
diamond makedb --in db.fasta --db ./database --taxonmap ../prot.accession2taxid --taxonnodes ../nodes.dmp --taxonnames ../names.dmp
mv database.dmnd ../db-wtax.dmnd

rm *.dmp readme.txt taxdump.tar.gz gc.prt

