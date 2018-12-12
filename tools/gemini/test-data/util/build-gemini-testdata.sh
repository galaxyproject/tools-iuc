cd "$(dirname "$0")"

export GEMINI_CONFIG=../test-cache
IN_PTH="$1"
OUT_PTH=$GEMINI_CONFIG/gemini/data
GENOMIC_REGION=3:187000000-187500000

# downsample all vcf and bed annotation files to the region of interest and reindex
for vcf in `ls $IN_PTH/*.gz | grep -v hprd_interaction_edges.gz -`
do
    python ./shrink_tabix.py $vcf -r $GENOMIC_REGION -o $OUT_PTH/`basename $vcf`
done

# downsample gene_table files to the region of interest
echo "$IN_PTH/summary_gene_table_v75 -> $OUT_PTH/summary_gene_table_v75"
python ./shrink_simple_tab.py $IN_PTH/summary_gene_table_v75 -r chr$GENOMIC_REGION -c 0 8 9 -n 1 -o $OUT_PTH/summary_gene_table_v75

echo "$IN_PTH/detailed_gene_table_v75 -> $OUT_PTH/detailed_gene_table_v75"
python ./shrink_simple_tab.py $IN_PTH/detailed_gene_table_v75 -r chr$GENOMIC_REGION -c 0 11 12 -n 1 -o $OUT_PTH/detailed_gene_table_v75

# filter kegg_pathway files to retain only records of the genes listed
# in the downsampled summary_gene_table
for kegg in `ls $IN_PTH/kegg_pathways_*`
do
    echo "$kegg -> $OUT_PTH/`basename $kegg`"
    cut -f2 $OUT_PTH/summary_gene_table_v75 | grep -Fv None | grep -Fwf - $kegg > $OUT_PTH/`basename $kegg`
done

# filter hprd_interaction file to retain only records of the genes listed
# in the downsampled summary_gene_table
echo "$IN_PTH/hprd_interaction_edges.gz -> $OUT_PTH/hprd_interaction_edges.gz"
bgzip -dc $IN_PTH/hprd_interaction_edges.gz > $OUT_PTH/hprd_interaction_edges
cut -f2 $OUT_PTH/summary_gene_table_v75 | grep -Fv None | grep -Ff - $OUT_PTH/hprd_interaction_edges | bgzip > $OUT_PTH/hprd_interaction_edges.gz
rm $OUT_PTH/hprd_interaction_edges

# filter cancer_gene_census file to retain only records of the genes listed
# in the downsampled summary_gene_table;
# TO DO: make the filter stricter by looking for matches only in the first
# column of the cancer_gene_census file (but the file is relatively small anyway)
echo "$IN_PTH/cancer_gene_census.20140120.tsv -> $OUT_PTH/cancer_gene_census.20140120.tsv"
cut -f2 $OUT_PTH/summary_gene_table_v75 | grep -Fv None | grep -Fwf - $IN_PTH/cancer_gene_census.20140120.tsv > $OUT_PTH/cancer_gene_census.20140120.tsv

# now use gemini load to build a test database
echo "Building gemini test database"
gemini --annotation-dir $OUT_PTH load -v ../gemini_load_input.vcf -t snpEff ../gemini_load_result1.db
gemini --annotation-dir $OUT_PTH load -v ../gemini_load_input.vcf -t snpEff --skip-gene-tables --no-load-genotypes ../gemini_load_result2.db
gemini --annotation-dir $OUT_PTH load -v ../gemini_load_input.vcf -t snpEff -p ../gemini_amend_input.ped ../gemini_load_result3.db

