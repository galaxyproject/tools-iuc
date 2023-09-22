# AlphaFold - for developers

The AlphaFold tool comprises several components:

- Wrapper `alphafold.xml` and associated `macro*.xml` files.
- Input FASTA validation `validate_fasta.py`
- Output file generation `outputs.py`
- Docker image - see `docker/` for building (hosted at hub.docker.com/neoformit/alphafold)
- AlphaFold mock for testing `fetch_test_data.sh`
- Test script for `outputs.py` - `tests/test_outputs.sh`
- Output data from different model presets for mocking - `test-data/*mer*_output/`
- HTML visualization output - `alphafold.html`


## Planemo testing

Add the following line to Planemo's Galaxy venv activation script (should be something like `~/.planemo/gx_venv_3/bin/activate`)

```sh
export PLANEMO_TEST=1
```

When you `planemo test` the wrapper should use the mock AlphaFold run, which copies AlphaFold outputs from `test-data/*mer*_output/` directories.


## Updating `outputs.py`

This script is a bit complex so it has it's own test script. This will test generation of outputs for the three model presets `monomer`, `monomer_ptm` and `multimer`:

```bash
# From the alphafold root dir
tests/test_outputs.sh

# Keep the test outputs to manually check them
tests/test_outputs.sh --keep
```
