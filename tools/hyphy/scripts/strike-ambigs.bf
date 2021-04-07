RequireVersion ("2.5.20");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");



filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment of coding sequences and replace any ambiguous codons with ---. Write results to a new file in FASTA format, and report changed sequences to stdout
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA"
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);

KeywordArgument ("code",      "Which genetic code should be used", "Universal");
KeywordArgument ("alignment", "An in-frame codon alignment in one of the formats supported by HyPhy");

filter.in =  alignments.PromptForGeneticCodeAndAlignment ("filter.dataset", "filter.input");

KeywordArgument     ("output", ".fasta for compressed data", None);
filter.out = io.PromptUserForFilePath(".fasta for filtered data");
fprintf (filter.out, CLEAR_FILE, KEEP_OPEN);

GetDataInfo (filter.site_patterns, filter.input);

filter.patter2site = {};


for (i,j,v; in; filter.site_patterns) {
    index = i+j;
    if (filter.patter2site / v == FALSE ) {
        filter.patter2site [v] = {};
    }  
    filter.patter2site [v] + index;
}

GET_DATA_INFO_RETURNS_ONLY_THE_INDEX = TRUE;
COUNT_GAPS_IN_FREQUENCIES = FALSE;
filter.unique_patterns = utility.Array1D (filter.input.site_freqs);

for (seq = 0; seq < filter.input.species; seq += 1) {
     io.ReportProgressBar ("filter","Processing sequence " + (1+seq));
     codons = {1, filter.input.sites};
     codons [0] = "";
     GetString (seq_name, filter.input, seq);
     GetDataInfo (seq_chars, filter.input, seq);

     filter.ambigs = 0;

     for (pattern = 0; pattern < filter.unique_patterns; pattern += 1) {
        GetDataInfo (pattern_info, filter.input, seq, pattern); 
        if (pattern_info >= 0) {
            codon_start = (filter.patter2site[pattern])[0] * 3;
            codon = seq_chars [codon_start][codon_start+2];
        } else {
            codon = "---";
            filter.ambigs += Abs (filter.patter2site [pattern])
        }
        for (c; in; filter.patter2site [pattern] ) {
            codons[c] = codon;
        }
     }
     if (filter.ambigs > 0) {
        fprintf (stdout, "\nStriking ", filter.ambigs, " codons that are incompletely resolved from " + seq_name + "\n");
     }
     fprintf (filter.out,">",seq_name,"\n",Join ("", codons), "\n");
}

fprintf (filter.out,CLOSE_FILE);


