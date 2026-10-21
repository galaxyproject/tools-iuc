# Test data for AMRFinderPlus

## Create `test-db`

The committed `test-db` is a reduced AMRFinderPlus database for the wrapper
tests. It is based on AMRFinderPlus database release `2026-05-15.1`
(`V4.2-2026-05-15.1` in `amrfinderplus_versioned.loc`) and keeps only the
content needed by the current tests.

1. Download the AMRFinderPlus database with `amrfinder_update`.
2. Keep in `AMRProt.fa` only sequences listed in `to_keep_for_test-db`.
3. Keep only the `Enterococcus_faecalis` row in `taxgroup.tsv`.
4. Keep only `Enterococcus` `AMR_DNA` files from the downloaded database.
5. Trim the retained protein and mutation metadata to the current test scope:
   keep only the header and `Enterococcus_faecalis` rows in
   `AMRProt-mutation.tsv`, keep only the header in
   `AMRProt-susceptible.tsv`, and keep `AMRProt-suppress.tsv` as an empty
   placeholder file.
6. Reduce `fam.tsv` to the hierarchy rows needed by the current test outputs
   and their parent rows.
7. Reduce `AMR.LIB` to one valid HMM model. Protein searches still require a
   syntactically valid HMM file, but the current tests do not depend on any HMM
   hits.
8. Run `amrfinder_index` on the reduced database directory:

    ```sh
    amrfinder_index test-db
    ```

9. ***Used AI agent to test different combinations of what can and cannot be removed or replaced***
Remove database files that are not used by the current wrapper tests:

    ```sh
    rm test-db/AMR.LIB.h3*
    rm test-db/AMR_CDS.fa*
    rm test-db/AMR_DNA-Enterococcus_faecalis.fa.njs
    rm test-db/AMR_DNA-Enterococcus_faecalis.fa.not
    rm test-db/AMR_DNA-Enterococcus_faecalis.fa.ntf
    rm test-db/AMR_DNA-Enterococcus_faecalis.fa.nto
    rm test-db/AMR_DNA-Enterococcus_faecium.fa*
    rm test-db/AMR_DNA-Enterococcus_faecium.tsv
    rm test-db/AMRProt.fa.pjs
    rm test-db/AMRProt.fa.ptf
    rm test-db/AMRProt.fa.pto
    ```

The final database keeps `AMR.LIB`, `AMRProt.fa` and the required `AMRProt.fa`
BLAST index files because protein searches still need them. The `AMR.LIB.h3*`
HMM index files are omitted because the tests pass with the reduced plain
`AMR.LIB`. `AMRProt-suppress.tsv` is kept as an empty file because combined
protein/nucleotide searches expect the file to exist. The
`Enterococcus_faecalis` mutation FASTA and required BLAST index files are kept
because the mutation-output tests exercise that path. Optional BLAST sidecar
files that are not needed by the tests are omitted. The `Enterococcus_faecium`
FASTA, BLAST index and mutation metadata files are omitted because no current
test uses faecium mutation search.
