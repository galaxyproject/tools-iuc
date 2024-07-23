## bigwig peak bed maker

### July 30 2024 for the VGP

This is a Galaxy tool, for building some of the [NIH MARBL T2T assembly polishing](https://github.com/marbl/training) tools as Galaxy workflows.

JBrowse2 2.12.3 update will include a plugin for optional colours to distinguish bed features, shown being tested in the screenshots below.

### Find and mark BigWig peaks to a bed file for display

In the spirit of DeepTools, but finding contiguous regions where the bigwig value is either above or below a given centile.
0.99 and 0.01 for example. These quantile cut point values are found and applied over each chromosome using some [cunning numpy code](http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html)

![image](https://github.com/fubar2/bigwig_peak_bed/assets/6016266/cdee3a2b-ae31-4282-b744-992c15fb49db)

![image](https://github.com/fubar2/bigwig_peak_bed/assets/6016266/59d1564b-0c34-42a3-b437-44332cf1b2f0)

Big differences between chromosomes 14,15,21,22 and Y in this "all contigs" view - explanations welcomed:

![image](https://github.com/fubar2/bigwig_peak_bed/assets/6016266/162bf681-2977-4eb8-8d6f-9dad5b3931f8)


[pybedtools](https://github.com/jackh726/bigtools) is used for the bigwig interface. Optionally allow
multiple bigwigs to be processed into a single bed - the bed features have the bigwig name in the label for viewing.

### Note on quantiles per chromosome rather than quantiles for the whole bigwig

It is just not feasible to hold all contigs in the entire decoded bigwig in RAM to estimate quantiles. It may be
better to sample across all chromosomes so as not to lose any systematic differences between them - the current method will hide those
differences unfortunately. Sampling might be possible. Looking at the actual quantile values across a couple of test bigwigs suggests that
there is not much variation between chromosomes but there's now a tabular report to check them for each input bigwig.

### Table reports

The optional table output report gives a crude histogram and the top/bottom 10 values to help 
understand what is likely to be informative. In this example, there are 26700 zero values so
using a lower cutoff quantile is likely to have a lot of them, although a large window requirement
will decease the overload...

Descriptive measures
bigwig  test
contig  chr10_PATERNAL
n       135711693
mean    12.178164
std     7.997467
min     0.000000
max     365.000000
qtop    364.00
qbot    noqlo
First/Last 10 value counts
Value   Count
0.00    26700
1.00    82900
2.00    261400
3.00    676993
4.00    1665500
5.00    3125700
6.00    5078000
7.00    7469000
8.00    10191700
9.00    12544600
355.00  100
356.00  100
357.00  300
358.00  100
360.00  500
361.00  300
362.00  200
363.00  600
364.00  900
365.00  700
Histogram of bigwig values
chr10_PATERNAL        18.25 | 127,047,593 | **************************************************************************
chr10_PATERNAL        36.50 |   7,510,000 | ****
chr10_PATERNAL        54.75 |     818,900 |
chr10_PATERNAL        73.00 |     117,200 |
chr10_PATERNAL        91.25 |      51,900 |
chr10_PATERNAL       109.50 |      44,200 |
chr10_PATERNAL       127.75 |      21,600 |
chr10_PATERNAL       146.00 |      17,900 |
chr10_PATERNAL       164.25 |      16,400 |
chr10_PATERNAL       182.50 |      18,600 |
chr10_PATERNAL       200.75 |       5,400 |
chr10_PATERNAL       219.00 |       6,600 |
chr10_PATERNAL       237.25 |       6,200 |
chr10_PATERNAL       255.50 |       3,900 |
chr10_PATERNAL       273.75 |       4,500 |
chr10_PATERNAL       292.00 |       7,100 |
chr10_PATERNAL       310.25 |       3,000 |
chr10_PATERNAL       328.50 |       2,700 |
chr10_PATERNAL       346.75 |       3,500 |
chr10_PATERNAL       365.00 |       4,500 |
chr10_PATERNAL ------------ |------------ |
chr10_PATERNAL           N= | 135,711,693 |
chr10_PATERNAL ------------ |------------ |


