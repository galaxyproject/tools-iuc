Introduction
============

`Kraken <http://ccb.jhu.edu/software/kraken/>`__ is a taxonomic sequence
classifier that assigns taxonomic labels to short DNA reads. It does
this by examining the :math:`k`-mers within a read and querying a
database with those :math:`k`-mers. This database contains a mapping of
every :math:`k`-mer in
`Kraken <http://ccb.jhu.edu/software/kraken/>`__'s genomic library to
the lowest common ancestor (LCA) in a taxonomic tree of all genomes that
contain that :math:`k`-mer. The set of LCA taxa that correspond to the
:math:`k`-mers in a read are then analyzed to create a single taxonomic
label for the read; this label can be any of the nodes in the taxonomic
tree. `Kraken <http://ccb.jhu.edu/software/kraken/>`__ is designed to be
rapid, sensitive, and highly precise. Our tests on various real and
simulated data have shown
`Kraken <http://ccb.jhu.edu/software/kraken/>`__ to have sensitivity
slightly lower than Megablast with precision being slightly higher. On a
set of simulated 100 bp reads,
`Kraken <http://ccb.jhu.edu/software/kraken/>`__ processed over 1.3
million reads per minute on a single core in normal operation, and over
4.1 million reads per minute in quick operation.

The latest released version of Kraken will be available at the `Kraken
website <http://ccb.jhu.edu/software/kraken/>`__, and the latest updates
to the Kraken source code are available at the `Kraken GitHub
repository <https://github.com/DerrickWood/kraken>`__.

If you use `Kraken <http://ccb.jhu.edu/software/kraken/>`__ in your
research, please cite the `Kraken
paper <http://genomebiology.com/2014/15/3/R46>`__. Thank you!

System Requirements
===================

Note: Users concerned about the disk or memory requirements should read
the paragraph about MiniKraken, below.

-  **Disk space**: Construction of Kraken's standard database will
   require at least 160 GB of disk space. Customized databases may
   require more or less space. Disk space used is linearly proportional
   to the number of distinct :math:`k`-mers; as of Feb. 2015, Kraken's
   default database contains just under 6 billion (6e9) distinct
   :math:`k`-mers.

   In addition, the disk used to store the database should be
   locally-attached storage. Storing the database on a network
   filesystem (NFS) partition can cause Kraken's operation to be very
   slow, or to be stopped completely. As NFS accesses are much slower
   than local disk accesses, both preloading and database building will
   be slowed by use of NFS.

-  **Memory**: To run efficiently, Kraken requires enough free memory to
   hold the database in RAM. While this can be accomplished using a
   ramdisk, Kraken supplies a utility for loading the database into RAM
   via the OS cache. The default database size is 75 GB (as of Feb.
   2015), and so you will need at least that much RAM if you want to
   build or run with the default database.

-  **Dependencies**: Kraken currently makes extensive use of Linux
   utilities such as sed, find, and wget. Many scripts are written using
   the Bash shell, and the main scripts are written using Perl. Core
   programs needed to build the database and run the classifier are
   written in C++, and need to be compiled using g++. Multithreading is
   handled using OpenMP. Downloads of NCBI data are performed by wget
   and in some cases, by rsync. Most Linux systems that have any sort of
   development package installed will have all of the above listed
   programs and libraries available.

   Finally, if you want to build your own database, you will need to
   install the
   `Jellyfish <http://www.cbcb.umd.edu/software/jellyfish/>`__
   :math:`k`-mer counter. Note that Kraken only supports use of
   Jellyfish version 1. Jellyfish version 2 is not yet compatible with
   Kraken.

-  **Network connectivity**: Kraken's standard database build and
   download commands expect unfettered FTP and rsync access to the NCBI
   FTP server. If you're working behind a proxy, you may need to set
   certain environment variables (such as ``ftp_proxy`` or
   ``RSYNC_PROXY``) in order to get these commands to work properly.

-  **MiniKraken**: To allow users with low-memory computing environments
   to use Kraken, we supply a reduced standard database that can be
   downloaded from the Kraken web site. When Kraken is run with a
   reduced database, we call it MiniKraken.

   The database we make available is only 4 GB in size, and should run
   well on computers with as little as 8 GB of RAM. Disk space required
   for this database is also only 4 GB.


