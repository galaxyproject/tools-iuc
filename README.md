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
