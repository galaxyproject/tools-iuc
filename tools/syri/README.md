# SyRI Galaxy wrapper

Galaxy tool wrapper for [SyRI](https://github.com/schneebergerlab/syri): Synteny and Rearrangement Identifier. Syri compares alignments between two chromosome-level assemblies and identifies synteny and structural rearrangements.

## Minimap2 Galaxy Wrapper Compatibility

* When using minimap2 for generating alignments for syri, the minimum required Galaxy wrapper version is `2.28+galaxy1`
* This version correctly adds the `--eqx` flag, which is required by SyRI to interpret alignment CIGAR strings

## Accepted Alignment Inputs

* For now alignments in `bam`, `sam` and `paf` are available as inputs
* Delta file from mummer is currently not accepted

## Sample name in output vcf

* Users can add a sample name to add to output vcf. This is checked by regex

```
  <param argument="--samplename" name="sample_name" type="text" label="Sample name for the output VCF file. (default: sample)" optional="true">
    <validator type="regex" message="Invalid characters in sample name">^[a-zA-Z0-9\-_]+$</validator>
  </param>
```

## Empty Map IDs Output — Not an Error

* The output **Map IDs file** lists corresponding chromosomes between the reference and query genomes
* SyRI only generates the output when the chromosome names differ between the two genomes.
* If the chromosome names are identical, the file will be **absent or empty** — this is expected and not an error.

## Exit codes

* Any non zero
* Regex pattern `- ERROR -` in stderr/stdout
