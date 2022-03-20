# SimText

A text mining framework for interactive analysis and visualization of similarities among biomedical entities.

## Brief overview of tools:

 - pubmed_by_queries: 

 For each search query, PMIDs or abstracts from PubMed are saved.

 - abstracts_by_pmids: 

 For all PMIDs in each row of a table the according abstracts are saved in additional columns.

 - text_to_wordmatrix: 

 The most frequent words of text from each row are extracted and united in one large binary matrix. 
 
 - pmids_to_pubtator_matrix: 

 For PMIDs of each row, scientific words are extracted using PubTator annotations and subsequently united in one large binary matrix. 

 - simtext_app: 

 Shiny app with word clouds, dimension reduction plot, dendrogram of hierarchical clustering and table with words and their frequency among the search queries.

## Set up user credentials on Galaxy

To enable users to set their credentials (NCBI API Key) for this tool,
make sure the file `config/user_preferences_extra_conf.yml` has the following section:

```
preferences:
    ncbi_account:
        description: NCBI account information
        inputs:
            - name: apikey
              label: NCBI API Key (available from "API Key Management" at https://www.ncbi.nlm.nih.gov/account/settings/)
              type: text
              required: False

```

## Requirements command-line version

 - R (version > 4.0.0)

## Installation command-line version

```
$ mkdir -p <path>/simtext
$ cd <path>/simtext
$ git clone https://github.com/dlal-group/simtext
```

## pubmed_by_queries

This tool uses a set of search queries to download a defined number of abstracts or PMIDs for each search query from PubMed. PubMed's search rules and syntax apply. Users can obtain an API key from the Settings page of their NCBI account (to create an account, visit http://www.ncbi.nlm.nih.gov/account/). If the tool is used as command-line tool the API key is passed as an argument. For usage in Galaxy the API key is added to the Galaxy user-preferences (User/ Preferences/ Manage Information).

Input:

Tab-delimited table with a list of search queries (biomedical entities of interest) in one column. The column header should start with "ID_" (e.g., "ID_gene" if search queries are genes). 

Usage:
```
$ Rscript pubmed_by_queries.R [-h] [-i INPUT] [-o OUTPUT] [-n NUMBER] [-a] [-k KEY] [--install_packages]
```

Optional arguments: 
```
 -h, --help                  show help message
 -i INPUT, --input INPUT     input file name. add path if file is not in working directory
 -o OUTPUT, --output OUTPUT  output file name [default "pubmed_by_queries_output"]
 -n NUMBER, --number NUMBER  number of PMIDs or abstracts to save per ID [default "5"]
 -a, --abstract              if abstracts instead of PMIDs should be retrieved use --abstracts 
 -k KEY, --key KEY           if NCBI API key is available, add it to speed up the download of PubMed data. For usage in Galaxy add the API key to the Galaxy user-preferences (User/ Preferences/ Manage Information).
 --install_packages          if you want to auto install missing required packages
```

Output: 

A table with additional columns containing PMIDs or abstracts from PubMed.

## abstracts_by_pmids

This tool retrieves abstracts for a matrix of PMIDs. The abstract text is saved in additional columns.

Input:

Tab-delimited table with rows representing biomedical entities and columns containing the corresponding PMIDs. The names of the PMID columns should start with “PMID_” (e.g., “PMID_1”, “PMID_2” etc.).

Usage:
```
$ Rscript abstracts_by_pmid.R [-h] [-i INPUT] [-o OUTPUT]
```

Optional arguments: 
```
 -h, --help                 show help message
 -i INPUT, --input INPUT    input file name. add path if file is not in working directory
 -o OUTPUT, --output OUTPUT output file name [default "abstracts_by_pmids_output"]
 --install_packages         if you want to auto install missing required packages
```

Output: 

A table with additional columns containing abstract texts.

## text_to_wordmatrix

The tool extracts for each row the most frequent words from the text in columns starting with "ABSTRACT" or "TEXT. The extracted words from each row are united in one large binary matrix, with 0= word not frequently occurring in text of that row and 1= word frequently present in text of that row.

Input: 

The output of ‘pubmed_by_queries’ or ‘abstracts_by_pmids’ tools, or a tab-delimited table with text in columns starting with "ABSTRACT" or "TEXT".

Usage:
```
$ Rscript text_to_wordmatrix.R [-h] [-i INPUT] [-o OUTPUT] [-n NUMBER] [-r] [-l] [-w] [-s] [-p]
```

Optional arguments: 
```
 -h, --help                    show help message
 -i INPUT, --input INPUT       input file name. add path if file is not in working directory
 -o OUTPUT, --output OUTPUT    output file name. [default "text_to_wordmatrix_output"]
 -n NUMBER, --number NUMBER    number of most frequent words that should be extracted per row [default "50"]
 -r, --remove_num              remove any numbers in text
 -l, --lower_case              by default all characters are translated to lower case. otherwise use -l
 -w, --remove_stopwords        by default a set of english stopwords (e.g., 'the' or 'not') are removed. otherwise use -w
 -s, --stemDoc                 apply Porter's stemming algorithm: collapsing words to a common root to aid comparison of vocabulary
 -p, --plurals                 by default words in plural and singular are merged to the singular form. otherwise use -p
 -- install_packages           if you want to auto install missing required packages
```

Output: 

A binary matrix in that each column represents one of the extracted words.

## pmids_to_pubtator_matrix

The tool uses all PMIDs per row and extracts "Gene", "Disease", "Mutation", "Chemical" and "Species" terms of the corresponding abstracts, using PubTator annotations. The user can choose from which categories terms should be extracted. The extracted terms are united in one large binary matrix, with 0= term not present in abstracts of that row and 1= term present in abstracts of that row. The user can decide if the scientific terms should be extracted and used as they are or if they should be grouped by their geneIDs/ meshIDs (several terms are often grouped into one ID). Also, by default all terms are extracted, otherwise the user can specify a number of most frequent words to extract per row.

Input: 

Output of 'abstracts_by_pmids' tool, or tab-delimited table with columns containing PMIDs. The names of the PMID columns should start with "PMID", e.g. "PMID_1", "PMID_2" etc.

Usage:
```
$ Rscript pmids_to_pubtator_matrix.R [-h] [-i INPUT] [-o OUTPUT] [-b BYID] [-n NUMBER][-c {Gene,Disease,Mutation,Chemical,Species} [{Gene,Disease,Mutation,Chemical,Species} ...]]
```
 
Optional arguments:
```
 -h, --help                    show help message
 -i INPUT, --input INPUT       input file name. add path if file is not in workind directory
 -o OUTPUT, --output OUTPUT    output file name. [default "pmids_to_pubtator_matrix_output"]
 -b, --byid                    if you want to find common gene IDs / mesh IDs instead of specific scientific terms.
 -n NUMBER, --number NUMBER    number of most frequent terms/IDs to extract. by default all terms/IDs are extracted.
 -c [...], --categories [...]  PubTator categories that should be considered [default "('Gene', 'Disease', 'Mutation','Chemical')"]
 -- install_packages           if you want to auto install missing required packages
```

Output: 

Binary matrix in that each column represents one of the extracted terms.

## simtext_app

The tool enables the exploration of data generated by ‘text_to_wordmatrix’ or ‘pmids_to_pubtator_matrix’ tools in a Shiny local instance. The following features can be generated: 1) word clouds for each initial search query, 2) dimension reduction and hierarchical clustering of binary matrices, and 3) tables with words and their frequency in the search queries.

Input:

1)	Input 1: 
Tab-delimited table with
	- A column with initial search queries starting with "ID_" (e.g., "ID_gene" if initial search queries were genes).
	- Column(s) with grouping factor(s) to compare pre-existing categories of the initial search queries with the grouping based on text. The column names should start with "GROUPING_". If the column name is "GROUPING_disorder", "disorder" will be shown as a grouping variable in the app.
2)	Input 2: 
The output of ‘text_to_wordmatrix’ or ‘pmids_to_pubtator_matrix’ tools, or a binary matrix.

Usage:
```
$ Rscript simtext_app.R [-h] [-i INPUT] [-m MATRIX] [-p PORT]
```

Optional arguments:
```
 -h,        --help             show help message
 -i INPUT,  --input INPUT      input file name. add path if file is not in working directory
 -m MATRIX, --matrix MATRIX    matrix file name. add path if file is not in working directory
 -p PORT,   --port PORT        specify port, otherwise randomly selected
 --host					specify host
 -- install_packages           if you want to auto install missing required packages
```

Output: 

SimText app
