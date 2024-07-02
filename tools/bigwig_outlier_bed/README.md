## bigwig peak bed maker

### July 30 2024 for the VGP

This code will soon become a Galaxy tool, for building some of the [NIH MARBL T2T assembly polishing](https://github.com/marbl/training) tools as Galaxy workflows.

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
