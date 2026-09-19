#!/usr/bin/env bash
declare -a listdb=("card" "resfinder" "megares" "argannot" "virulencefinder" "plasmidfinder")

for i in "${listdb[@]}"
do
   ariba getref $i $i
   ariba prepareref -f ${i}.fa -m ${i}.tsv ${i}.db
   rm *.{fa,tsv,log}
done
