A suite of 14 tools for the **[Montreal Forced Aligner (MFA)](https://montreal-forced-aligner.readthedocs.io/en/latest/)**, a state-of-the-art command-line utility for aligning speech and text using Kaldi.

### Tools Included
This suite wraps the core MFA functionality, split into logical categories:

* **Alignment:** `mfa_align`
* **Validation:** `mfa_validate`, `mfa_validate_dictionary`, `mfa_find_oovs`
* **Transcription (ASR):** `mfa_transcribe` (Kaldi)
   - `mfa_transcribe_whisper` (Whisper) is excluded as we have Whisper as a separate tool on Galaxy
* **Preparation Utilities:** `mfa_g2p`, `mfa_tokenize`, `mfa_remap_alignments`, `mfa_model_add_words`
  - `mfa_remap_dictionary` is excluded as I could not perform a functional test using a small sample locally.
* **Training & Adaptation:** `mfa_train`, `mfa_adapt`, `mfa_train_dictionary`, `mfa_train_g2p`, `mfa_train_lm`, `mfa_train_tokenizer`
  - `mfa_train_ivector` is excluded for now as I could not wrap it properly
* **Diarization:** `mfa_diarize`

### Technical Implementation & Scope

* **Data Tables:** This implementation uses standard 4-column Data Tables (`value`, `dbkey`, `name`, `path`) to manage Pretrained Models (Acoustic, Dictionary, G2P, IVector, LM, Tokenizer). *As this is my first implementation of Galaxy Data Tables, I would appreciate a specific review of the `.loc` file generation logic and table configuration.*
* **Scope:** This suite covers the most important commands. Some subcommands were omitted either because they strictly require GPU resources (not available on all standard Galaxy nodes) or due to complexity in generating valid unit test data for edge cases.
* **Domain Context:** Please note that speech processing is different from my background. While the tools pass technical validation and linting, feedback from a domain expert regarding parameter defaults and descriptions would be highly valued. However, this is only possible if we have it on usegalaxy.eu for our specific users to test them.

### Instructions for Reviewers & Admins

To test this tool locally or deploy it, you **should** set up the reference data first. The tools rely on pre-downloaded models to avoid massive downloads during job execution.

**1. Setup Directory & Download Models**
Edit and run `download_mfa_models.sh`. You likely need to adjust the target directory variable to match your server's storage:

In download_mfa_models.sh
TARGET_DIR="/data/db/mfa"  # <--- Change this to your preferred path
Run: `bash download_mfa_models.sh [GITHUB_TOKEN]`

**2. Generate Data Tables Edit and run generate_loc_files.py to create the Galaxy .loc files based on the models downloaded in step 1.**
In generate_loc_files.py
MFA_MODELS_ROOT = "/data/db/mfa/pretrained_models" # <--- Must match the path from step 1
Run `python3 generate_loc_files.py`

**3. Generate Test Data To run planemo test, generate the lightweight dummy datasets (sine waves and synthetic text).**
Run: python3 generate_test_data.py

### Testing Strategy

- Synthetic Data: I created a Python script (generate_test_data.py) included in the repo to generate valid dummy data without copyright issues.
- Optimization: Tests use small dictionaries (4 words) and robust noise-injected audio to ensure MFA's mathematical validation (LDA/MLLT) passes quickly without requiring massive datasets.