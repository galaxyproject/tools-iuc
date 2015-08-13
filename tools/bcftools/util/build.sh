#!/bin/bash

# Convert
for x in help/*.txt; do 
    python util/convert.py $x > bcftools_$(basename $x .txt | sed 's/help_//g').xml; 
done;
# Auto-macros
python util/macroify.py --macros macros.xml bcftools_*.xml
# Undo changes to macros file because don't know how to add conditionally.
git checkout -- macros.xml

# Prettify
for i in *.xml; do cat $i | xmllint --pretty 1 - > $i.pretty; mv $i.pretty $i; python util/reformat.py $i; done;

