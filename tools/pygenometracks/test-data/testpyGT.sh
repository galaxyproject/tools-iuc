. <(planemo conda_env pyGenomeTracks.xml)
pgt --tracks test-data/test1.ini --region chrX:3000000-3500000 -o test-data/master_TADs_plot.png
pgt --tracks test-data/test2.ini --region chrX:3000000-3500000 -o test-data/bigwig_multiple.png
pgt --tracks test-data/test3.ini --region chrX:3000000-3500000 -o test-data/master_TADs_BW_plot.png
pgt --tracks test-data/test4.ini --region X:2700000-3100000 -o test-data/test_alpha.png
pgt --tracks test-data/test5.ini --region X:3000000-3300000 -o test-data/test_gtf_bed4.png
pgt --tracks test-data/test6.ini --region X:2760000-2802000 -o test-data/test_narrowPeak.png
pgt --tracks test-data/test7.ini --region chrX:3300000-3500000 -o test-data/test_gtf_flybase_param.png
pgt --tracks test-data/test8.ini --region chrX:3300000-3500000 -o test-data/test_ucsc_param.png
pgt --tracks test-data/test9.ini --region X:3130000-3140000 -o test-data/test_arrowhead_zoom.png
pgt --tracks test-data/test10.ini --region X:3340000-3380000 -o test-data/test_middle_triangle.png
pgt --tracks test-data/test11.ini --region chrX:3250000-3400000 -o test-data/test_TADs_bdgm.png

conda_env_deactivate
