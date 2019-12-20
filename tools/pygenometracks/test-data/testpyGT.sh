. <(planemo conda_env pyGenomeTracks.xml)
pgt --tracks test-data/test1.ini --region chrX:3000000-3500000 -o test-data/master_TADs_plot.png
pgt --tracks test-data/test2.ini --region chrX:3000000-3500000 -o test-data/bigwig_multiple.png
pgt --tracks test-data/test3.ini --region chrX:3000000-3500000 -o test-data/master_TADs_BW_plot.png
pgt --tracks test-data/test4.ini --region X:2700000-3100000 -o test-data/test_alpha.png
pgt --tracks test-data/test5.ini --region X:3000000-3300000 -o test-data/test_gtf_bed4.png
pgt --tracks test-data/test6.ini --region X:2760000-2802000 -o test-data/test_narrowPeak.png

