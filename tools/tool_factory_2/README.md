toolfactory_2
=============

This is an upgrade to the tool factory but with added parameters 
(optionally editable in the generated tool form - otherwise fixed) and 
multiple input files.

Any number of  parameters up to the limit of your patience with repeat groups
These are optionally editable by the user - names cannot be changed so 
no overwriting $JAVA_HOME_DIR or else permanently fixed and not editable at run time.

Any number of input files can be passed to your script, but of course it
has to deal with them. Both path and metadata name are supplied either in the environment 
(bash/sh) or as command line parameters (python,perl,rscript) that need to be parsed and
dealt with in the script.

