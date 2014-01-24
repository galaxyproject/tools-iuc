These are galaxy tools for SnpEff a variant annotation and effect prediction tool by Pablo Cingolani. 
It annotates and predicts the effects of variants on genes (such as amino acid changes).
( http://snpeff.sourceforge.net/ )

This repository contains a tool_dependencies.xml file that will attempt to automatically install SnpEff and SnpSift.   

This will use the default location for genome reference downloads from the snpEff.config:
data_dir = ~/snpEff/data/
You can manually edit the installed snpEff.config and change the location, or you can create a symbolic link to the desired data location from ~/snpEff.

The genome reference options used by the tools:
    "SnpEff"  snpEff.xml
    "SnpEff Download" snpEff_download.xml
are taken from: tool-data/snpeffect_genomedb.loc

There are 2 datamanagers to download and install prebuilt SnpEff Genome databases:
  data_manager_snpeff_databases - generates a list of available SnpEff genome databases into the tool-data/snpeff_databases.loc 
  data_manager_snpeff_download - downloads a SnpEff genome database selected from: tool-data/snpeff_databases.loc and adds entries to snpeff_genomedb.loc,snpeff_regulationdb.loc,snpeff_annotations.loc 

SnpEff citation:
"A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672 [PubMed - in process]

SnpSift citation:
"Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift", Cingolani, P., et. al., Frontiers in Genetics, 3, 2012.

-Wrapper authors:
    -- Jim Johnson

