Galaxy wrapper for Augustus
===========================

This wrapper is copyright 2012-2013 by Björn Grüning.

This is a wrapper for the command line tool of Augustus_.

.. _augustus: http://bioinf.uni-greifswald.de/augustus/

AUGUSTUS is a program that predicts genes in eukaryotic genomic sequences.

Oliver Keller, Martin Kollmar, Mario Stanke, Stephan Waack (2011)
A novel hybrid gene prediction method employing protein multiple sequence alignments
Bioinformatics, doi: 10.1093/bioinformatics/btr010

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008)
Using native and syntenically mapped cDNA alignments to improve de novo gene finding
Bioinformatics, doi: 10.1093/bioinformatics/btn013

Mario Stanke and Stephan Waack (2003)
Gene Prediction with a Hidden-Markov Model and a new Intron Submodel.
Bioinformatics, Vol. 19, Suppl. 2, pages ii215-ii225


History
=======

- v0.1: Initial public release
- v0.2: Added tool_dependencies.xml file and update the augustus version (thanks to James Johnson)
- v0.3: upgrade to augustus 2.7, added new organisms and new parameters, output additional sequence files
- v0.3.1: added parallelism and changed the output parameters from boolean to a select box
- v3.2.3: updated to augustus 3.2.3, use conda dependency, migrated to IUC
