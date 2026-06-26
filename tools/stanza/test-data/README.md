# Test Data for Stanza NLP

## Running Tests Locally

The tests require Stanza language models to be available locally. The easiest way is to copy from an existing Stanza installation:

```bash
# If you have Stanza models cached locally:
cp -r /Users/$USER/Library/Caches/stanza/1.11.0/resources test-data/stanza_models

# OR download models using Python:
python3 -c "
import stanza
stanza.download('en', model_dir='test-data/stanza_models', package='default_fast')
"
```

Then run tests from the tool directory:

```bash
planemo test --docker
```

## Required Models

The tool uses the `default_fast` package which includes nocharlm models for:
- tokenize
- pos (part-of-speech tagging)  
- ner (named entity recognition)
- depparse (dependency parsing)
- sentiment

## Data Table Configuration

The `test-data/stanza_models.loc` file points to `${__HERE__}/stanza_models` which resolves to the local test-data directory during testing.

## Model Size

The complete model directory is ~1.2GB and is excluded from the repository via `.gitignore`.

## Docker Integration  

With the `--docker` flag, Galaxy automatically mounts the `test-data/` directory into the container, making the models accessible to Stanza at runtime.