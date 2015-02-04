# Comparative Genomics Circos tools for use in Galaxy

# Spec

- [ ] 1v1 genome plots showing regions with high levels of synteny (color coded)
    - [ ] write Mauve/progressiveMauve library to handle processing of data and identifying syntenic regions
    - [ ] write data processing tool that takes this and prepares it for use in Circos
    - [ ] (challenge) write a post-processor that can annotate the image with information on hover with useful facts like geneinformation for the genes on the plots, synteny information for the synteny links.
- [ ] attempt plots with more genomes exploring intergenic relationships, exploring if they’re an effective tool
- [ ] explore [Hive plots](http://www.hiveplot.net/) (another thing by the circos people) to see if we can effectively represent the required information in it.
