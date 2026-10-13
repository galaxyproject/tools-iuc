# Beacon2 Galaxy wrapper drafts

Files included:
- `macros.xml`
- `csv_to_bff.xml`
- `genomicVariations_vcf.xml`
- `genomicVariations_postprocessing.xml`
- `individuals_to_cohorts_csv.xml`

## Design notes

These wrappers are designed around the current Python CLI exposed by Beacon2 RI Tools v2:
- `csv_to_bff.py --output --datasetId --input --entry_type`
- `genomicVariations_vcf.py --output --datasetId --refGen --caseLevelData --numRows --verbosity --json --input --alleleCounts --alleleFrequency`
- `genomicVariations_postprocessing.py --input --datasetId`
- `individuals_to_cohorts_csv.py --input --output --datasetId --cohortId --cohortName --cohortType`

The wrappers generate a temporary `conf/conf.py` at runtime and copy the installed package tree into the job directory before execution. This is meant to cope with the upstream scripts relying on package-relative resources such as `files/headers/*.txt`, `files/deref_schemas/*.json`, and `conf/conf.py`.

## Important caveat

I could validate that the XML files are well-formed, but I did **not** run `planemo test` here. You will still want to adapt test-data paths and, if needed, the exact conda package name/version used by your Bioconda recipe.
