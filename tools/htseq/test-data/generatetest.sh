#python ../htseqsams2mx.py -g rn4_chr20_100k.gtf -o test.xls --samf "'rn4chr20test1.bam','col1'" --samf "'rn4chr20test2.bam','col2'"
python ../htseqsams2mx.py -g rn4_chr20_100k.gtf -o htseqsams2mx_test1_out.xls --samf "'rn4chr20test1.bam','test1','bam','rn4chr20test1.bai'" --samf "'rn4chr20test2.bam','test2','bam','rn4chr20test2.bai'"

