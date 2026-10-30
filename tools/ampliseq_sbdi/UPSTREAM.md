# Upstream export formatters

The two R scripts, normalized with IUC styler, are from nf-core/ampliseq 2.11.0, commit `0473e157ac9a7b1d36845ff9f8fa7ac843d3b234`, `bin/`. Original author: Daniel Lundin. MIT copyright and permission notices are retained in each R source file.

The wrappers normalize the first ASV/taxonomy/prediction column header to ASV_ID. Metadata must retain the original ID column. Annotation provenance fields are supplied by the user; these exporters do not classify sequences. date_identified records execution date and is not a scientific measurement.

Galaxy metadata are generated separately from the actual output headers by galaxy_metadata.R. The formatter calculations are unchanged. The declared dplyr dependency is needed for across(); container selection is left to Galaxy dependency resolution.
