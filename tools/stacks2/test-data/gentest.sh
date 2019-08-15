#!/usr/bin/env bash

mkdir stacks_outputs

denovo_map.pl --samples demultiplexed --popmap denovo_map/popmap_cstacks.tsv -o stacks_outputs --paired  && 
gunzip -c stacks_outputs/catalog.calls > stacks_outputs/catalog.calls.vcf
rm stacks_outputs/catalog.calls

mv stacks_outputs/PopA_0{1,2}.{tags,snps,alleles}.tsv ustacks/
mv stacks_outputs/catalog.{tags,snps,alleles}.tsv cstacks/
mv stacks_outputs/PopA_0{1,2}.matches.tsv sstacks/
mv stacks_outputs/PopA_0{1,2}.matches.bam tsv2bam/
mv stacks_outputs/tsv2bam.log tsv2bam/
mv stacks_outputs/catalog.{calls.vcf,fa.gz} gstacks/
mv stacks_outputs/gstacks.log* gstacks/
mv stacks_outputs/populations.* populations/
mv stacks_outputs/denovo_map.log denovo_map/

rmdir stacks_outputs
