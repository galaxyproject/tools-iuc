# Test data for AMRFinderPlus

## Create `test-db` 

1. Download ARMFinderPlus database from NCBI using the data-manager Python script after modifying the `main` function to keep only `amrfinderplus_download.download_amrfinderplus_db()`
2. Remove `AMR_DNA` files that are not `Enterococcus`
3. Remove lines not related to `Enterococcus` in `taxgroup.tab`
4. Keep in `AMR_CDS`, `AMRProt`, and `ReferenceGeneCatalog.txt` only sequences listed in `to_keep_from_test-db`
5. Keep in `amr_targets.fa` only sequences listed in `amr_targets_to_keep_from_test-db`
6. Copy the `AMR.LIB` from a previous version
7. Run the data-manager Python script after modifying `main` function like this:

    ```
    amrfinderplus_download.download_amrfinderplus_db()
    amrfinderplus_download.amrfinderplus_db_path = f'{amrfinderplus_download._output_dir}/{amrfinderplus_download._db_name}'
    amrfinderplus_download.make_hmm_profile()
    amrfinderplus_download.make_blastdb()
    ```

6. Move the `amrfinderplus-db` folder content to `test-db`