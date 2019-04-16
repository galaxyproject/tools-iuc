Introduction
============
Kraken is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. Kraken examines the $k$-mers within a query sequence and uses the information within those $k$-mers to query a database. That database maps $k$-mers to the lowest common ancestor (LCA) of all genomes known to contain a given $k$-mer.

The first version of Kraken used a large indexed and sorted list of $k$-mer/LCA pairs as its database. While fast, the large memory requirements posed some problems for users, and so Kraken 2 was created to provide a solution to those problems.

Kraken 2 differs from Kraken 1 in several important ways:

1. Only minimizers of the $k$-mers in the query sequences are used as database queries. Similarly, only minimizers of the $k$-mers in the reference sequences in the database's genomic library are stored in the database. We will also refer to the minimizers as $\ell$-mers, where $\ell \leq k$. All $k$-mers are considered to have the same LCA as their minimizer's database LCA value.
2. Kraken 2 uses a compact hash table that is a probabilistic data structure. This means that occasionally, database queries will fail by either returning the wrong LCA, or by not resulting in a search failure when a queried minimizer was never actually stored in the database. By incurring the risk of these false positives in the data structure, Kraken 2 is able to achieve faster speeds and lower memory requirements. Users should be aware that database false positive errors occur in less than 1% of queries, and can be compensated for by use of confidence scoring thresholds.
3. Kraken 2 has the ability to build a database from amino acid sequences and perform a translated search of the query sequences against that database.
4. Kraken 2 utilizes spaced seeds in the storage and querying of minimizers to improve classification accuracy.
5. Kraken 2 provides support for "special" databases that are not based on NCBI's taxonomy. These are currently limited to three popular 16S databases.

Because Kraken 2 only stores minimizers in its hash table, and $k$ can be much larger than $\ell$, only a small percentage of the possible $\ell$-mers in a genomic library are actually deposited in the database. This creates a situation similar to the Kraken 1 "MiniKraken" databases; however, preliminary testing has shown the accuracy of a reduced Kraken 2 database to be quite similar to the full-sized Kraken 2 database, while Kraken 1's MiniKraken databases often resulted in a substantial loss of per-read sensitivity.

The Kraken 2 paper is currently under preparation. Until it is released, please cite the original Kraken paper if you use Kraken 2 in your research. Thank you!
Page: https://ccb.jhu.edu/software/kraken2/

System Requirements
===================
- Disk space: Construction of a Kraken 2 standard database requires approximately 100 GB of disk space. A test on 01 Jan 2018 of the default installation showed 42 GB of disk space was used to store the genomic library files, 26 GB was used to store the taxonomy information from NCBI, and 29 GB was used to store the Kraken 2 compact hash table.

- Like in Kraken 1, we strongly suggest against using NFS storage to store the Kraken 2 database if at all possible.

- Memory: To run efficiently, Kraken 2 requires enough free memory to hold the database (primarily the hash table) in RAM. While this can be accomplished with a ramdisk, Kraken 2 will by default load the database into process-local RAM; the --memory-mapping switch to kraken2 will avoid doing so. The default database size is 29 GB (as of Jan. 2018), and you will need slightly more than that in RAM if you want to build the default database.

- Dependencies: Kraken 2 currently makes extensive use of Linux utilities such as sed, find, and wget. Many scripts are written using the Bash shell, and the main scripts are written using Perl. Core programs needed to build the database and run the classifier are written in C++11, and need to be compiled using a somewhat recent version of g++ that will support C++11. Multithreading is handled using OpenMP. Downloads of NCBI data are performed by wget and rsync. Most Linux systems will have all of the above listed programs and development libraries available either by default or via package download.

- Unlike Kraken 1, Kraken 2 does not use an external $k$-mer counter. However, by default, Kraken 2 will attempt to use the dustmasker or segmasker programs provided as part of NCBI's BLAST suite to mask low-complexity regions (see [Masking of Low-complexity Sequences]).

- MacOS NOTE: MacOS and other non-Linux operating systems are not explicitly supported by the developers, and MacOS users should refer to the Kraken-users group for support in installing the appropriate utilities to allow for full operation of Kraken 2. We will attempt to use MacOS-compliant code when possible, but development and testing time is at a premium and we cannot guarantee that Kraken 2 will install and work to its full potential on a default installation of MacOS.

- In particular, we note that the default MacOS X installation of GCC does not have support for OpenMP. Without OpenMP, Kraken 2 is limited to single-threaded operation, resulting in slower build and classification runtimes.

- Network connectivity: Kraken 2's standard database build and download commands expect unfettered FTP and rsync access to the NCBI FTP server. If you're working behind a proxy, you may need to set certain environment variables (such as ftp_proxy or RSYNC_PROXY) in order to get these commands to work properly.

- Kraken 2's scripts default to using rsync for most downloads; however, you may find that your network situation prevents use of rsync. In such cases, you can try the --use-ftp option to kraken2-build to force the downloads to occur via FTP.

- MiniKraken: At present, users with low-memory computing environments can replicate the "MiniKraken" functionality of Kraken 1 in two ways: first, by increasing the value of $k$ with respect to $\ell$ (using the --kmer-len and --minimizer-len options to kraken2-build); and secondly, through downsampling of minimizers (from both the database and query sequences) using a hash function. This second option is performed if the --max-db-size option to kraken2-build is used; however, the two options are not mutually exclusive. In a difference from Kraken 1, Kraken 2 does not require building a full database and then shrinking it to obtain a reduced database.


