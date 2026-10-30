# Upstream export formatters

The two R scripts, normalized with IUC styler, are from nf-core/ampliseq 2.11.0, commit `0473e157ac9a7b1d36845ff9f8fa7ac843d3b234`, `bin/`. Original author: Daniel Lundin. MIT license retained in LICENSE.upstream.

The wrappers normalize the first ASV/taxonomy/prediction column header to ASV_ID. Metadata must retain the original ID column. Annotation provenance fields are supplied by the user; these exporters do not classify sequences. date_identified records execution date and is not a scientific measurement.
