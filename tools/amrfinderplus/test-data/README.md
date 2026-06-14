# Test data for AMRFinderPlus

## Create `test-db`

The committed `test-db` is a reduced AMRFinderPlus database for the wrapper
tests. It is based on AMRFinderPlus database release `2026-05-15.1`
(`V4.2-2026-05-15.1` in `amrfinderplus_versioned.loc`) and keeps only the
content needed by the current tests.

1. Download the AMRFinderPlus database with `amrfinder_update`.
2. Keep in `AMRProt.fa` only sequences listed in `to_keep_for_test-db`.
3. Keep only `Enterococcus` rows in `taxgroup.tsv`.
4. Keep only `Enterococcus` `AMR_DNA` files from the downloaded database.
5. Trim the retained protein and mutation metadata to the targets listed in
   `amr_targets_to_keep_for_test-db`, where applicable.
6. Copy the reduced `AMR.LIB` from a previous test database version.
7. Run `amrfinder_index` on the reduced database directory:

    ```sh
    amrfinder_index test-db
    ```

8. Remove database files that are not used by the current wrapper tests:

    ```sh
    rm test-db/AMR.LIB.h3*
    rm test-db/AMR_CDS.fa*
    rm test-db/AMR_DNA-Enterococcus_faecium.fa*
    rm test-db/AMR_DNA-Enterococcus_faecium.tsv
    ```

The final database keeps `AMR.LIB`, `AMRProt.fa` and the `AMRProt.fa` BLAST
index files because protein searches still need them. The `AMR.LIB.h3*` HMM
index files are omitted because the tests pass with the reduced plain
`AMR.LIB`. The `Enterococcus_faecalis` mutation FASTA and BLAST index files are
kept because the mutation-output tests exercise that path. The
`Enterococcus_faecium` FASTA, BLAST index and mutation metadata files are
omitted because no current test uses faecium mutation search.
