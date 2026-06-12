# Test data for AMRFinderPlus

## Create `test-db` 

1. Download the AMRFinderPlus database with `amrfinder_update`.
2. Remove `AMR_DNA` files that are not `Enterococcus`.
3. Remove lines not related to `Enterococcus` in `taxgroup.tsv`.
4. Keep in `AMR_CDS.fa` and `AMRProt.fa` only sequences listed in `to_keep_for_test-db`.
5. Copy the reduced `AMR.LIB` from a previous version.
6. Run `amrfinder_index` on the reduced database directory:

    ```
    amrfinder_index test-db
    ```
