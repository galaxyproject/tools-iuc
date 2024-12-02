ResFinder Database documentation
=============

The ResFinder database is a curated database of acquired resistance genes.

## Content of the repository
1. <AMR class>.fsa - DNA sequence of the genes from a specified AMR class
2. phenotypes.txt - Translation table for genotype to phenotype prediction
3. phenotype_panels.txt - For selected species, lists the relevant phenotypes
4. INSTALL.py - Script for indexing the database with KMA
5. CHECK-entries.sh - Script examining whether there is consistency between the databases and the phenotypes.txt
6. antibiotic_classes.txt - Lists the relevant antibiotic classes
7. config - Configuration fie that describes the content of the databases
8. notes.txt - Lists additional notes of the entries in the databases
9. history.txt - history file of ResFinderDB commits dating back to 2018/01/01

## Installation
Download the ResFinder database by cloning the database:
```bash
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
```
The database can be used with BLAST as-is.

If you want to use the database with the stand-alone ResFinder tool, and wishes
to use the mapping based method (available from ResFinder version 4.0.0), the
database needs to be indexed with KMA.

### Installing KMA (optional):

If you are running the stand-alone ResFinder in docker, you may be able to skip
installing KMA, and just rely on the temporary KMA installation done by the
Docker script.

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
cd path/to/resfinder_db
python3 INSTALL.py /path/to/kma_index non_interactive
```
The path to kma_index can be omitted if it exists in PATH or if the script
should attempt to do an automatic temporary installation of KMA.

#### c) Index database manually (not recommended)
It is possible to index the databases manually, but is generally not recommended
as it is more prone to error. If you choose to do so, be aware of the naming of
the indexed files.

This is an example of how to index the ResFinder database files:
```bash
# Go to the database directory
cd path/to/resfinder_db
# create indexing directory
mkdir kma_indexing
# Index files using kma_index
kma_index -i fusidicacid.fsa -o kma_indexing/fusidicacid
kma_index -i phenicol.fsa -o kma_indexing/phenicol
kma_index -i glycopeptide.fsa -o kma_indexing/glycopeptide
kma_index -i trimethoprim.fsa -o kma_indexing/trimethoprim
kma_index -i oxazolidinone.fsa -o kma_indexing/oxazolidinone
kma_index -i tetracycline.fsa -o kma_indexing/tetracycline
kma_index -i quinolone.fsa -o kma_indexing/quinolone
kma_index -i nitroimidazole.fsa -o kma_indexing/nitroimidazole
kma_index -i fosfomycin.fsa -o kma_indexing/fosfomycin
kma_index -i aminoglycoside.fsa -o kma_indexing/aminoglycoside
kma_index -i macrolide.fsa -o kma_indexing/macrolide
kma_index -i sulphonamide.fsa -o kma_indexing/sulphonamide
kma_index -i rifampicin.fsa -o kma_indexing/rifampicin
kma_index -i colistin.fsa -o kma_indexing/colistin
kma_index -i beta-lactam.fsa -o kma_indexing/beta-lactam
```

### Documentation

The documentation available as of the date of this release can be found at:
https://bitbucket.org/genomicepidemiology/resfinder_db/overview.


Citation
=======

When using the method please cite:

Identification of acquired antimicrobial resistance genes.
Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup 
FM, Larsen MV.
J Antimicrob Chemother. 2012 Jul 10.
PMID: 22782487         doi: 10.1093/jac/dks261

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
