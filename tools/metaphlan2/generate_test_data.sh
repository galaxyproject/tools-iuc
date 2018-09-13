#/usr/bin/env bash

metaphlan_hclust_heatmap.py \
    --in test-data/merged_community_profile.tabular \
    --out test-data/heatmap.png \
    -m 'average' \
    -d 'braycurtis' \
    -f 'correlation' \
    --minv '0' \
    --tax_lev 'a' \
    --sdend_h '0.1' \
    --fdend_w '0.1' \
    --cm_h '0.03' \
    --font_size '7' \
    --clust_line_w '1' \
    --perc '90' \
    -c 'jet'

metaphlan_hclust_heatmap.py \
    --in test-data/merged_community_profile.tabular \
    --out test-data/heatmap.pdf \
    -m 'ward' \
    -d 'euclidean' \
    -f 'euclidean' \
    --minv '0' \
    --tax_lev 'a' \
    --sdend_h '0.1' \
    --fdend_w '0.1' \
    --cm_h '0.03' \
    --font_size '7' \
    --clust_line_w '1' \
    --perc '90' \
    -c 'pink'

metaphlan_hclust_heatmap.py \
    --in test-data/merged_community_profile.tabular \
    --out test-data/heatmap.svg \
    -m 'complete' \
    -d 'hamming' \
    -f 'matching' \
    --minv '0' \
    --tax_lev 'a' \
    --sdend_h '0.1' \
    --fdend_w '0.1' \
    --cm_h '0.03' \
    --font_size '7' \
    --clust_line_w '1' \
    --perc '90' \
    -c 'pink'