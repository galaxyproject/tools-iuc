binary_file=ex.bed
target_markers=6000
target_cols=$(( 6 + ( $target_markers * 2 )))

dec_prefix=out
subs_prefix=sub
new_prefix="new."${target_markers}

##conda activate plink
## Extract binary file into new map, ped, and fam files
plink --bfile $(basename $binary_file .bed) --recode --out ${dec_prefix}

## Subset ped and map
cat ${dec_prefix}.ped | cut -d' ' -f 1-${target_cols} > ${subs_prefix}.ped
head ${dec_prefix}.map -n ${target_markers} > ${subs_prefix}.map

## Recode into new
plink --file ${subs_prefix} --out ${new_prefix}
echo ${new_prefix}*
##conda deactivate
