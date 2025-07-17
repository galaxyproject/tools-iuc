# A proposal for tool panel restructuring and tool review

> This is a draft for IUC consideration

The amount of tools is overwhelming and we were not doing a very consistent job in maintaining tool sections, keeping subsections, sorting tools, and ensuring that each tool has a truly comprehensive help section. This is OK since the tool set has been (and is) undergoing rapid expansion. But we do need to get it under control as it interferes with fundamental usability of Galaxy.

## Table of contents

- Tool panel naming conventions
- Help section style
- Tool categories

## Tool panel naming conventions

The tool panel is narrow and long names do not work well. Tool names need to be to the point, so all unnecessary words, articles, punctuation signs must be avoided. All this information should be provided in the `<help>` section instead.

Authors of new tools should specify which category their tool should go to. The tools should ideally be tagged to facilitate search and alternative rendering of the tool panel for field-specific instances or individual users.

### Proposed guidelines

- Tool names should be as short and as expressive as possible 
- The bold portion of the tool name (defined by the `name` attribute within the `<tool>` tag) should be short 
- No parentheses are allowed in tool names
- Only the first word is capitalized

### Examples

| Current name | Suggested edit |
|--------------|---------------|
|**Select first** lines from a dataset | **Select first** lines |
|**Sort** data in ascending and descending order | **Sort** data|
|**Samtools fastx** extract FASTA or FASTQ from alignment files | SAMTools should have a dedicated category within the tool panel (see **Tool categories** below). The suggested name then would be: "**fastx**: extract FASTA or FASTQ"|

## Help section style

The Help section of tools vary widely in content. However, just like we cannot accept tools without tests, we should not be accepting tools with rudimentary tool sections. 

### Proposed guidelines

These guidelines were partially inspired by the structure of `man` pages:

- Each tool should include `Synopsis` section with a brief, one line phrase ending with a period. It prescribes what the tool does and is not a tool description.
- `Synopsis` should be followed by `Description` section
- If this is a homegrown tool, the `Description` section can be free form
- If this is an external package (e.g., SAMTools), the `Description` section should be a copy of help section for this tool adopted for Galaxy. This means that it should only describe CLI flags that are present within `<inputs>` section.
- When appropriate, the help section should include images preferably in [SVG](https://en.wikipedia.org/wiki/Scalable_Vector_Graphics) format.
- For complex tools the `<help>` section should include examples.

### Examples

Before:

![](https://i.imgur.com/oV4uY5m.png)

After:

![](https://i.imgur.com/G6rxKq3.png)

## Tool Categories

Tool section needs to be completely reworked. The existing section dates back to 2012/2013 and no longer reflects the breadth of Galaxy's toolset. 

### Proposed tool panel structure

Here top level sections (e.g., `Get Data`) are labels such as, for example, **GENERAL TEXT TOOLS** on usegalaxy.org. The next level are expandable tool categories containing actual tools. One open question is how to handle tool from different disciplines (e.g., genomics vs. climate science). Within each section categories and tools must be alphabetically listed.

- Get data
    - EBI
    - NCBI SRA
    - NCBI Datasets
    - UCSC Table Browsers
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
- Statistics and ML
    -  IWTomics
    -  ML
- Format conversion
    - FASTA/FASTQ
    - BED, GFF
- Sequence data 
    - FASTA manipulation
    - SeqTk
    - FASTQ manipulation
    - QC
    - EMBOSS
    - Fast5
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
- RNA-seq
    - Mapping
    - Transcript assembly
    - Quantification
    - SingleCell
    - RSeQC
- ChIP-seq
    -  Peak calling
    -  Motif analysis
- Immunology
    - Presto
- Proteomics
    - MS data analysis
    - Docking modeling
- Computational chemistry
    - Cheminformatics
    - Protein-ligand docking
    - Molecular dynamics
- Evolution
    - Phylogenetic tree reconstruction
    - Selection analysis
    - PlantTribes
- Climate science
    - [A Climate expert needs to set categories here]
