PointFinder Database documentation
=============

The PointFinder database is a curated database of resistance causing
chromosomal point mutations.

## Content of the repository
_TODO_

## Installation
Clone the database
```bash
git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git
```
The database can be used with BLAST as-is.

If you want to use the database with the stand-alone ResFinder tool, and wishes
to use the mapping based method (available from ResFinder version 4.0.0), the
database needs to be indexed.

### Installing KMA (optional):

If you are running the stand-alone ResFinder in docker, you may be able to skip
installing KMA, and just rely on the temporary KMA installation done by the
INSTALL script (or if you are just lazy and don't want to type the path to kma).

If you are not running ResFinder stand-alone in docker, you will need to
install KMA, if the mapping based method is needed (recommended).

#### Download and install KMA
```bash
# Go to the directory in which you want KMA installed
cd /some/path
# Clone KMA
git clone https://bitbucket.org/genomicepidemiology/kma.git
# Go to kma directory and compile code
cd kma && make
```

### Indexing with *INSTALL.py*
If you have KMA installed you either need to have the kma_index in your PATH or
you need to provide the path to kma_index to INSTALL.py

#### a) Run INSTALL.py in interactive mode
```bash
# Go to the database directory
cd path/to/resfinder_db
python3 INSTALL.py
```
If kma_index was found in your path a lot of indexing information will be
printed to your terminal, and will end with the word "done".

If kma_index wasn't found you will recieve the following output:
```bash
KMA index program, kma_index, does not exist or is not executable
Please input path to executable kma_index program or choose one of the options below:
	1. Install KMA using make, index db, then remove KMA.
	2. Exit
```
You can now write the path to kma_index and finish with <enter> or you can
enter "1" or "2" and finish with <enter>.

If "1" is chosen, the script will attempt to install kma in your systems
default temporary location. If the installation is successful it will proceed
to index your database, when finished it will delete the kma installation again.

#### b) Run INSTALL.py in non_interactive mode
```bash
# Go to the database directory
cd path/to/pointfinder_db
python3 INSTALL.py /path/to/kma_index non_interactive
```
The path to kma_index can be omitted if it exists in PATH or if the script
should attempt to do an automatic temporary installation of KMA.

#### c) Index database manually (not recommended)
It is possible to index the databases manually, but is generally not recommended
as it is more prone to error. If you choose to do so, be aware of the naming of
the indexed files.

This is an example of how to index the PointFinder database files:
```bash
# Go to the database directory
cd path/to/pointfinder_db
# create indexing directory
mkdir kma_indexing
# Index files using kma_index
kma_index -i db_pointfinder/campylobacter/*.fsa -o db_pointfinder/campylobacter/campylobacter
kma_index -i db_pointfinder/escherichia_coli/*.fsa -o db_pointfinder/escherichia_coli/escherichia_coli
kma_index -i db_pointfinder/enterococcus_faecalis/*.fsa -o db_pointfinder/enterococcus_faecalis/enterococcus_faecalis
kma_index -i db_pointfinder/enterococcus_faecium/*.fsa -o db_pointfinder/enterococcus_faecium/enterococcus_faecium
kma_index -i db_pointfinder/neisseria_gonorrhoeae/*.fsa -o db_pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae
kma_index -i db_pointfinder/salmonella/*.fsa -o db_pointfinder/salmonella/salmonella
kma_index -i db_pointfinder/mycobacterium_tuberculosis/*.fsa -o db_pointfinder/mycobacterium_tuberculosis/mycobacterium_tuberculosis
```

## PointFinder database format

Each species that the PointFinder database covers has its own folder in which the database for the corresponding species resides.
Four types of files exists within a PointFinder database:

1. One or more FASTA files ending with the extension ".fsa". These FASTA files contains the reference (wild type) sequence of a specific region, mutations are found with respect to these sequences. The fasta header in these files must contain just the gene name as given in the resistens-overview.txt file.
2. A file called "genes.txt". It is used to describe which of the FASTA sequences that should be employed when using the database.
3. One or no file called "RNA_genes.txt". Defines which of the FASTA files are RNA genes.
4. A file called "resistens-overview.txt". Defines the resistance causing mutations and phenotypes. The File is described in details below.

### resistens-overview.txt

The file is a text file in tab separated format. The first line starts with a #, followed by the headers for the table. Indels are always described at the end. If any indels are described, a line consisting of only "# Indels", should precede the first indel entry (row).

|     Header   | Explanation                                                                                                       |
| -------------|-------------------------------------------------------------------------------------------------------------------|
| Gene_ID      | Gene ID as written in the genes.txt file                                                                          |
| Gene_name    | Name of gene or region                                                                                            |
| Codon_pos    | Nucleotide position in the FASTA file where the mutation starts                                                   |
| Ref_nuc      | The reference sequence at the mutation, "-" if mutation is an indel                                               |
| Ref_codon    | One letter aa or nucleotide describing the reference sequence. Can also be "del" or "ins" if mutation is an indel |
| Res_codon    | Comma separated list of nucleotides or amino acids (1-letter code)                                                |
| Resistance   | Comma separated list of antibiotics                                                                               |
| PMID         | Comma separated list of pubmed IDs describing the mutation                                                        |
| Mechanism    | Description of the resistance mechanism                                                                           |
| Notes        | Text with other information                                                                                       |
| Required_mut | Other mutations needed in order to gain resistance (see below for more details)                                   |

### Required_mut format

There are several layers that needs to be addressed in this field. Starting at the bottom, if a mutation is only dependent on a single other mutation it is written like so:
```
<Gene_ID>_<Ref_codon><Codon_pos><muts>
```
Where Gene_ID, Ref_codon, and Codon_pos are described in the above table. "muts" are either a single mutation written in 1-letter code or a number of possible mutations separated by a dot. Like so:
```
<Gene_ID>_<Ref_codon><Codon_pos><mut1>.<mut2>.<mut3>...<mutN>
```
Example of a required mutation in gyrA at position 83, changing an S to either L, W, A, or V:
```
gyrA_S83L.W.A.V
```
If there are several required mutations, that all need to be present, they are separated by commas, like so:
```
M = <Gene_ID>_<Ref_codon><Codon_pos><muts>
M,M,M
```
Example:
```
pmrA_S39I,pmrA_R81S
```
If there are several groups of mutations that each can confer resistance with the mutation in question, but is independent of each other, the groups are separated by semicolons, like so:
```
M1,M1,M1;M2,M2,M2
```
Note that a group can consist of just one required mutation.
Example of a required mutation in gyrA either at position 83 or 87:
```
gyrA_S83Y.F.A;gyrA_D87N.G.Y.K
```

## Documentation

The documentation available as of the date of this release can be found at
https://bitbucket.org/genomicepidemiology/pointfinder_db/overview.


Citation
=======

When using the method please cite:

Not yet published


License
=======

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
