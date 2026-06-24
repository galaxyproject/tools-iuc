# Galaxy wrapper for spaCy

This wrapper provides easy access to spaCy NLP pipelines with multi-language support. Users can select language models, annotators, and output formats. See the [spaCy website](https://spacy.io) for more information.

## Features

- **75+ Languages**: Support for a wide range of languages with trained models
- **Multiple Model Sizes**: Choose between small (fast), medium (balanced), large (accurate), or transformer-based (most accurate) models
- **Flexible Annotation**: Tokenization, POS tagging, NER, dependency parsing, and more
- **Multiple Output Formats**: JSON, CoNLL, CoNLL-U, and human-readable text

## Supported Languages

Models are available for 25+ languages including:

- English (en)
- Spanish (es)
- German (de)
- French (fr)
- Chinese (zh)
- Japanese (ja)
- Korean (ko)
- Portuguese (pt)
- Italian (it)
- Dutch (nl)
- Greek (el)
- Polish (pl)
- Norwegian Bokmål (nb)
- Lithuanian (lt)
- Danish (da)
- Swedish (sv)
- Romanian (ro)
- Catalan (ca)
- Finnish (fi)
- Croatian (hr)
- Macedonian (mk)
- Russian (ru)
- Ukrainian (uk)

**Note:** Not all annotators are available for all languages. See the [spaCy models documentation](https://spacy.io/usage/models) for language-specific capabilities.

## Installation

### For Galaxy Administrators

This tool requires language models to be installed via the data manager:

1. **Configure Data Tables**: Copy `tool_data_table_conf.xml.sample` to your Galaxy configuration
2. **Create Location File**: Copy `tool-data/spacy_models.loc.sample` to `tool-data/spacy_models.loc`
3. **Install Language Models**: Use the "spaCy Language Models" data manager to download and register models
   - Navigate to Admin → Local Data in Galaxy
   - Run the data manager for each model you want to support
   - Models are installed via pip using `python -m spacy download`

### Manual Installation

Alternatively, you can manually install models and add entries to `spacy_models.loc`:

```bash
# Install a spaCy model
python -m spacy download en_core_web_sm

# Add entry to tool-data/spacy_models.loc (tab-separated)
en_core_web_sm	English (small)	en	sm	en_core_web_sm
```

## Model Sizes

spaCy models come in different sizes with trade-offs between speed and accuracy:

- **sm (small)**: 10-20MB, fastest processing, good for high-volume tasks
- **md (medium)**: 40-50MB, balanced performance, includes word vectors
- **lg (large)**: 500-800MB, highest accuracy for traditional models, large word vectors
- **trf (transformer)**: 400-500MB, transformer-based, state-of-the-art accuracy

## Annotation Types

### Tokenization and sentence segmentation
Splits text into tokens and identifies sentence boundaries.

### Part of speech and lemmas
Includes tokenization plus part-of-speech (POS) tagging and lemmatization (reducing words to their base form).

### Named entity recognition (NER)
Includes tokenization, POS, and lemmas, plus identification of named entities such as PERSON, ORG, GPE, DATE, etc.

### Dependency parsing
Includes all previous annotators plus syntactic dependency parsing to identify grammatical relationships between words.

## Output Formats

- **JSON**: Comprehensive structured output with all annotations. Best for programmatic access and further processing.
- **CoNLL**: Tab-separated format suitable for dependency parsing tasks. Standard format used in NLP research.
- **CoNLL-U**: Universal Dependencies format with morphological features. More detailed than CoNLL.
- **Text**: Human-readable text output with statistics and annotations. Good for inspection and debugging.

## Comparison with Stanford CoreNLP

### Advantages of spaCy:
- **Faster processing**: Optimized for production use
- **More languages**: 75+ languages vs 8 in CoreNLP
- **Easier installation**: Pure Python, installable via pip
- **Modern architecture**: Built for modern NLP workflows
- **Better documentation**: Extensive examples and tutorials

### When to use CoreNLP:
- Need specific Stanford NLP features like coreference resolution
- Working with existing CoreNLP pipelines
- Require Java-based tools integration

## LICENSE

This tool wrapper is provided under the MIT License.

## Citation

When using spaCy, please cite:

```bibtex
@software{spacy,
  author = {Honnibal, Matthew and Montani, Ines and Van Landeghem, Sofie and Boyd, Adriane},
  title = {{spaCy: Industrial-strength Natural Language Processing in Python}},
  year = 2020,
  publisher = {Zenodo},
  doi = {10.5281/zenodo.1212303},
  url = {https://spacy.io}
}
```
