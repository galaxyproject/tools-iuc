# GFF3 validator 0.1.0 wrapper

Executable wrapper, macros and fixtures imported unchanged from
FAIR-bioHeaders/gff3-validator commit
`fda9715f2e776ce8453ba4d1fe9fd070951d8e88`, directory `galaxy/`.
The Tool Shed repository URL points to the existing IUC submission branch
containing this copy, so repository metadata matches the submitted directory.
The canonical upstream project remains the homepage URL.
The included LICENSE is copied exactly from that same source commit: applicable
current contributions retain MPL-2.0 and earlier MIT grants remain intact.
This directory does not change the license of other tools-iuc material or the
published validator package.

The wrapper targets the published Bioconda `gff3-validator=0.1.0`. It offers
plain/gzip GFF3, optional plain/gzip genome FASTA and translation-table biology
checks, header modes, HTML and optional JSON/SARIF reports. Invalid-file exit1
is a successful report-producing job; incomplete-run exit2 fails the job.

The release has no `--profile` option; main-only profiles and the later lone-CR
reader correction are not offered or claimed by this wrapper. Reports retain
explicit skipped/planned checks, including incomplete SO coverage. A model or
header result is not biological proof or Sequence Ontology endorsement.

Verification before submission: strict Planemo0.75.32 lint; eight complete
Galaxy25.0/Python3.12 tests passed with the noneditable published0.1.0 CLI
(0failures,0errors,0skips). Source executable wrapper, macros and all fixtures remain byte-identical to
the reviewed upstream source; only the repository-location metadata is adapted
for this IUC submission.
