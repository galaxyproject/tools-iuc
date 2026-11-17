## microsatellites to bed features

 **Convert short repetitive sequences to bed features**

 Microsatellites are usually defined as repeated short DNA patterns in an unbroken sequence.
 A microsatellite pattern or *motif* can be any combination nucleotides, typically from 1 to 6nt in length.
 
 This tool allows microsatellite and related features to be selected from a fasta sequence input file, and output into a single bed track, suitable for viewing in a genome browser such as JBrowse2.

 All motifs of selected lengths can be reported as individual features in the output bed file, or specific motifs can be provided and all
 others will be ignored. In all cases, a minimum required number of repeats can be specified. For example, requiring 2 or more repeats of the trimer *ACG* will report 
 every sequence of *ACGACG* or *ACGACGACG* or *ACGACGACGACG* and so on, as individual bed features.  Similarly, requiring 3 repeats of any trimer will 
 report every distinct 3 nucleotide pattern, including *ACGACGACG* as well as every other unique 3 nucleotide pattern with 3 sequential repeats or more such, as "CTCCTCCTC*.

 For other output formats, the pytrf native command line *findstr* can be used to produce a gff, csv or tsv output containing all exact short tandem repeats, as 
 described at the end of https://pytrf.readthedocs.io/en/latest

 A fasta file must be supplied for processing. A built in genome can be selected, or a fasta file of any kind can be selected from the current history. Note that all 
 symbols are treated as valid nucleotides by pytrf, so extraneous characters such as *-* or *N* in the input fasta may appear as unexpected bed features. Lower case fasta symbols will be converted
 to uppercase, to prevent them being reported as distinct motifs.


 **Filter motifs by length**
 
 The default tool form setting is to select all dimer motif patterns. 
 
 Additional motif lengths from 1 to 6nt can be selected in the multiple-select drop-down list. All features will be returned in a single bed file. For each selected motif length, 
 the minimum number of repeats required for reporting can be adjusted. **Tandem repeats** are defined as at least 2 of any pattern. This tool allows singleton motifs to be reported,
 so is not restricted to short tandem repeats (STR)

 **Filter motifs by pattern**

 This option allows a motif pattern to be specified as a text string such as *CG* or *ATC*. Multiple motifs can be specified as a comma separated string such as *CG,ATC*.
 All features will be returned as a single bed file.

 The minimum number of repeats for all motifs can be set to match specific requirements.

 For example, technical sequencing read bias may be influenced by the density of specific dimers, whether they are repeated or not
 such as in https://github.com/arangrhie/T2T-Polish/tree/master/pattern
 
 **Run pytrf findstr to create a csv, tsv or gff format output with all perfect STR**

This selection runs the pytrf *findstr* option to create gff/csv/tsv outputs as described at the end of https://pytrf.readthedocs.io/en/latest/. 

Quote below

*A Tandem repeat (TR) in genomic sequence is a set of adjacent short DNA sequence repeated consecutively. The core sequence or repeat unit is generally called a motif. 
According to the motif length, tandem repeats can be classified as microsatellites and minisatellites. Microsatellites are also known as simple sequence repeats (SSRs) 
or short tandem repeats (STRs) with motif length of 1-6 bp. Minisatellites are also sometimes referred to as variable number of tandem repeats (VNTRs) has longer motif length than microsatellites.
Pytrf is a lightweight Python C extension for identification of tandem repeats. The pytrf enables to fastly identify both exact or perfect SSRs.
It also can find generic tandem repeats with any size of motif, such as with maximum motif length of 100 bp. Additionally, it has capability of finding approximate or imperfect tandem repeats*

 