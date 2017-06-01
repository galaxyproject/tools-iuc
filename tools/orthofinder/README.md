Orthofinder only_groups for Galaxy

Version 1.0 - June 2017

  - This tool allows to run the first part of Orthofinder (Emms, D.M. and Kelly, S., 2015).
    -> runs blastp + MCL to infer orthogroups
  - Conda dependencies for Orthofinder 1.1.4 and python 2.7.

*******************************************************************************

This readme is currently used for implementation details ; it will be re-writed correctly afterwoods.

June 7th 2017 WORK IN PROGRESS 

  - Tool splitted in two part
  - 1st part is finished : infers orthogroups (blastp+MCL)
    implement arguments : -f, -b, -og, -t, -a
  - 2d part is about to start

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

  - Outputs with multiples subfiles and subfolders. Not easy to implement.
