#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

if (($#))
then
	genome=$1
	cut -f -4 < g1.G.bed > group1.bed
	bedGraphToBigWig group1.bed <(fetchChromSizes "$genome") group1.bw
	cut -f -4 < g2.G.bed > group2.bed
	bedGraphToBigWig group2.bed <(fetchChromSizes "$genome") group2.bw
fi
