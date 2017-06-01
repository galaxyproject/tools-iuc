Changelog

This readme is currently used for implementation details ; it will be re-writed correctly afterwoods.

01/06/2017 - Version 1.0 - WORK IN PROGRESS 

  - Started to split the tool in two part
  - 1st part is finished : it goes up to the -og option : It infers orthogroups (blastp+MCL)
    Available options : -f, -b, -op, -og, -t, -a
  - The 2d part is about to start

*******************************************************************************

26/05/2017 - Version 1.0 - WORK IN PROGRESS
  
  - Conda dependency for Orthofinder 1.1.4
  - The wrapper already includes the following tool options : -f, -b, -fg, -op, -og, -I, -s, -t, -a
  - Missing options : 
     - -ft (quite hard to implement : it needs inputs with subdirectories and files from all the previous steps), 
     - -x
  - No readme for the following options : -os, -oa, -ot, -M
  - Tests are made with the Example Dataset of Orthofinder (4 species of Mycoplasma, quick enought to make many tests).
  - The wrapper is NOT (yet?) implemented to deal with incompatible options. The user has to know a bit what he is doing.

  # BUGs and/or Issues : #

  - About the outputs : according to the selected option, the directory for the output files regularly differs (the tool is written this way in the 1.1.4 version, maybe it has been improved int the recent versions): sometimes, a "Results_MonthDay" directory is created, sometimes not. Acording to this, I had to add multiple if and elif statements in the command, in order to have the outputs files always at the same place for the galaxy outputs.

  - For now, the input files are not send in the tool as a dataset collection : the history quickly becomes messy when using options like -b or -fg (many files required).
  - Issue about dataset collections for outputs : Some files are missing, and sub-directories are ignored ! Once the bug resolved, I will try to improve that and write the output in order to have only one collection with all the results.
    - Example of missing file : galaxy does not seem to make the difference between "Orthogroups.txt" and "Orthogroups.csv" : only .txt file is in the collection.
    - Already Fixed (not sure) : at some point, the working directory output was incomplete (some mising files) when running orthofinder from the very beginning ("orthofinder -f .").
