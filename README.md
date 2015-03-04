# Comparative Genomics Circos tools for use in Galaxy

The CPT often needs to do 1v1 and 1vN genome comparisons (NvN are more rare and we have some other tools for that volume of information). These tools will support this sort of plotting to enable biologists and end-users to more easily visualize their data.

## Spec

- [ ] 1v1 genome plots showing regions with high levels of synteny (color coded)
    - [ ] write data processing tool that takes this and prepares it for use in Circos
        - progressiveMauve is available at the command line on the cpt server
        - Use it to align some test sequences, see how it behaves.
        - based on the output of progressiveMauve, write a tool to parse those files and produce a Circos Plot (once you get here bug me for what to do about Circos plots)
    - [ ] (challenge) write a post-processor that can annotate the image with information on hover with useful facts like geneinformation for the genes on the plots, synteny information for the synteny links.
- [ ] attempt plots with more genomes exploring intergenic relationships, exploring if they’re an effective tool
- [ ] explore [Hive plots](http://www.hiveplot.net/) (another thing by the circos people) to see if we can effectively represent the required information in it.

Circos Notes

	Configuration File
		The configuration file is the central element of a Circos graphic.
		Configuration files can be used to direct the creation of PNG or SVG images.
		Configuration files are parsed with the Config::General module.
		Files can include blocks for ideograms, links, etc.
		There are many ways to specify colors: see http://circos.ca/documentation/tutorials/configuration/configuration_files/ for details.
		Separate data files are required for chromosomes, tracks, links, and highlights.
	Karyotype File
		The following is an example of a Karyotype data file.
			chr - hs1 1 0 249250621 chr1
			chr - hs2 2 0 243199373 chr2
			chr - hs3 3 0 198022430 chr3
		The columns of this file represent genome name, label, start position, end position, and color.	
		Most Circos stuff is chromosome-based -- can links be drawn between tracks?
	Links
		Below is an example of a "links" block to be passed to Circos.
			<links>
			<link>
			#Each <link> block contains a different link data file.
			#Per line in a link file, different link-specific options can be specified
			#Link specific options should be tab-delimited from the data and comma-separated from each other at the end of each link in the data.
			file          = data/5/segdup.txt
			color         = black_a5
			radius        = 0.95r
			bezier_radius = 0.1r
			thickness     = 1
			#Number of records read from the data file can be capped with the parameter record_limit.
			</link>
			</links>
		For the size of the link to reflect the sizes of the matches on the two "chromosomes", the ribbon setting must be set to "yes" in the conf file.

