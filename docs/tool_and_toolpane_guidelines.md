# A proposal for tool panel restructuring and tool review

> This is a draft for IUC consideration

The amount of tools is overwhelming and we were not doing a very consistent job in maintaining tool sections, keeping subsections, sorting tools, and ensuring that each tool has a truly comprehensive help section. This is OK since the tool set has been (and is) undergoing rapid expansion. But we do need to get it under control as it interferes with fundamental usability of Galaxy.

## Table of contents

- Tool panel naming conventions
- Help section style
- Tool categories

## Tool panel naming conventions

The tool panel is narrow and long names do not work well. Tool names need to be to the point, so all unnecessary words, articles, punctuation signs must be avoided. All this information should be provided in the `<help>` section instead.

Author of every new tool should specify which category her/his tool should go to. The tools should ideally be tagged to facilitate search and alternative rendering of tool panel for field-specific instances or individual users.

### Proposed guidelines

- Tool names should be as short and as expressive as possible 
- The bold portion of the tool name (defined by `name` attribute within `<tool>` tag) should be short 
- No parenthesis are allowed in tool names
- Only the first word is capitalized

### Examples

| Current name | Suggested edit |
|--------------|---------------|
|**Select first** lines from a dataset | **Select first** lines |
|**Sort** data in ascending and descending order | **Sort** data|
|**Samtools fastx** extract FASTA or FASTQ from alignment files | SAMTools should have a dedicated category within the tool panel (see **Tool categories** below). The suggested name then would be: "**fastx**: extract FASTA or FASTQ"|

## Help section style

Help section of the tools vary widely in content. However, just like we cannot accept tools without tests we should not be accepting tools with rudimentary tool sections. 

### Proposed guidelines

These guidelines partially inspired by `man` pages structure:

- Each tool should include `Synopsis` section with a brief, one line phrase ending with a period. It prescribes what the tool does and is not a tool description.
- `Synopsis` should be followed by `Description` section
- If this is a homegrown tool, the `Description` section can be free form
- If this is an external package (e.g., SAMTools), the `Description` section should be a copy of help section for this tool adopted for Galaxy. This means that it should only describe CLI flags that are present within `<inputs>` section.
- When appropriate help section should include images preferably in `svg` format.
- For complex tools `<help>` section should include examples.

### Examples

Before:

![](https://i.imgur.com/oV4uY5m.png)

After:

![](https://i.imgur.com/G6rxKq3.png)

## Tool Categories

Tool section needs to be completely reworked. The existing section dates back to 2012/2013 and no longer reflects the breadth of Galaxy's toolset. 

### Proposed tool panel structure

Here top level sections (e.g., `Get Data`) are labels such as, for example, **GENERAL TEXT TOOLS** on usegalaxy.org. The next level are expandable tool categories containing actual tools. One open question is how to handle tool from different disciplines (e.g., genomics vs. climate science).

- Get data
    - NCBI SRA
    - NCBI Datasets
    - UCSC Table Browsers
    - EBI
- Interactive tools
    - Notebook environments
    - Other IEs
- Dataset collections
    - Collection operations
- Text and dataframes
    - Filtering [filters of all sorts]
    - Find/Replace [sed ...]
    - Column manipulation [cut1, cut2, paste, merge...]
    - Sorting [sort]
    - Join, Subtract, Group
    - Misc text tools
    - Statistical operations
    - Datamash
    - IWTomics
- Format conversion
    - FASTA/FASTQ
    - BED, GFF
- FASTA/FASTQ
    - FASTA manipulation
    - SeqTk
    - FASTQ manipulation
    - QC
    - EMBOSS
- Genomic intervals
    - BEDTools
- Read mapping
    - Short read data
    - Long read data
- Mapped data
    - SAMTools
    - Picard
    - BEDTools [can we list same category twice?]
    - DeepTools
    - QualiMap
    - UMI-tools
- Assembly
    - Viral and bacterial
    - Animals and plants
    - Annotation tools
    - VGP tools
- Metagenomics
    - Kraken
    - dada2
    - Mothur
    - VSearch
    - Diamond
    - Vegan
    - Diamond
    - CAT
    - Misc 
- Variation data
    - VCFTools
    - BCFTools
    - SNPEff
    - Diploid callers
    - Mixed/Haploid callers
    - MiModD
    - DuNovo
- RNAseq
    - Mapping
    - Transcript assembly
    - Quantification
    - SingleCell
    - RSeQC
-  ChIP-seq
    -  Peak calling
    -  Motif analysis
- Immunology
    - Presto
- Protein analysis
    - MS data analysis
    - Docking modeling
- Machine learning
    - [A ML expert will need to set categories here]
- Computational chemistry
    - Chemical tool box
- Evolution
    - Phylogenetic tree reconstruction
    - Selection analysis
    - PlantTribes
- Climate science
    - [A Climate expert needs to set categories here]
