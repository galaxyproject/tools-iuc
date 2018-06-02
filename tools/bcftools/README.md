# bcftools (v1.3)

Copied from branch bcftools1.2:

This aims to be a "faithful" rendering of the bcftool suite. I.e. options are
presented essentially as closely to the command line version as is useful.

This may not appeal to all, if you'd like to see smaller and more dedicated
tools (e.g. "intersect", "union" and "complement" being separate tools instead
of all of them included in the "isec" tool,) please feel free to file an issue.

Updated for bcftools v1.3

This was extended from the bcftools1.2 branch then greatly hand edited to
group params and manage param innteractions.

In the macros.xml there are macros and tokens to handle file input and output.
These use the datatypes currently available in galaxy: Vcf and Bcf
The macros take care of bgzip and indexing of inputs.

The convert command was split into 2 tools, "convert to vcf" and "convert from vcf"

## TODO:

- stats needs a matplotlib tool dependency  and pdflatex for generating a pdf of plots
- cnv needs a matplotlib tool dependency for generating images, then a means to consolidate those.
- cnv needs an input.vcf for testing, runs with bcftools cnv -s "HG00101" -o 'HG00101/' -p 5 mpileup.vcf
- roh needs a more useful input.vcf for testing
- plugin color chrs
- plugin frameshifts

## Status

The wrappers were automatically generated in bulk. That doesn't get them 100%
of the way there (e.g. meaningful test cases), so the rest of the process is a
bit slower.

- [x] annotate
- [x] call
- [ ] cnv (needs real test data, needs plotting)
- [x] concat
- [x] consensus
- [x] convert from vcf
- [x] convert to vcf
- [x] filter
- [x] gtcheck
- [x] isec
- [x] merge
- [x] norm
- [x] query
- [x] query list samples
- [x] reheader
- [x] roh
- [x] stats (needs plotting)
- [x] view
- [ ] +color chrs
- [x] +counts
- [x] +dosage
- [x] +fill an ac
- [x] +fill tags
- [x] +fixploidy
- [ ] +frameshifts
- [x] +impute info
- [x] +mendelian
- [x] +missing2ref
- [x] +setgt
- [x] +tag2tag
- [x] +vcf2sex
