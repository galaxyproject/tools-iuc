# Galaxy Data Manager for Stanza Language Models

This Galaxy data manager downloads and installs Stanza language models for use with the Stanza NLP annotation tool, supporting a curated selection of Stanza's 80+ languages with neural models trained on Universal Dependencies.

## Features

- **Curated language selection**: Install one or more of ~50 commonly used languages
- **Selectable model package**: Choose `default_fast`, `default`, or `default_accurate` depending on the annotators and accuracy you need
- **Constituency parsing support**: The `default` and `default_accurate` packages add constituency parsing (absent from `default_fast`)
- **Multiple packages per language**: The same language can be installed under more than one package; each becomes a separate selectable entry in the tool
- **Data table integration**: Registers models in Galaxy's `stanza_models` data table for the Stanza NLP tool

## How It Works

This data manager:
1. **Downloads with Stanza**: Uses the `stanza.download()` API to fetch models for the selected package directly from Stanford's HuggingFace repository.
2. **Installs per language**: Each language is stored as a self-contained `stanza_resources` directory (a `resources.json` plus a per-language model subdirectory).
3. **Registers models**: Writes an entry to the `stanza_models` data table so the Stanza NLP tool can load it with `stanza.Pipeline(dir=models_path, package=<package>)`.
4. **Managed storage**: Models are downloaded into the job's writable output directory, then moved by Galaxy into the managed data directory (see `data_manager_conf.xml`).
5. **Version**: Installs models compatible with Stanza 1.12.0.

## Model Packages

| Package | Contents | Constituency? | Approx. English size |
|---------|----------|---------------|----------------------|
| `default_fast` | Memory-efficient nocharlm models: tokenization, POS, lemma, depparse, NER, sentiment (where available) | No | ~130–200 MB |
| `default` | charlm-based models; adds constituency parsing | Yes | ~300 MB |
| `default_accurate` | Largest, most accurate (electra-based) models; adds constituency parsing | Yes | Largest |

`default_fast` is recommended for most deployments. Install `default` or `default_accurate` only for languages where you need constituency parsing or maximum accuracy, as they are substantially larger.

## Supported Languages

The data manager offers a curated selection of ~50 of Stanza's 80+ languages, including:

### European Languages
- **Western**: English, German, French, Spanish, Italian, Portuguese, Dutch
- **Nordic**: Swedish, Danish, Norwegian (Bokmål/Nynorsk), Finnish
- **Slavic**: Russian, Ukrainian, Polish, Czech, Slovak, Croatian, Serbian, Bulgarian
- **Other**: Greek, Hungarian, Romanian, Estonian, Latvian, Lithuanian, Catalan, Slovenian

### Asian Languages
- **East Asian**: Chinese (Simplified/Traditional), Japanese, Korean
- **South Asian**: Hindi, Tamil, Telugu, Marathi, Urdu
- **Southeast Asian**: Vietnamese, Thai, Indonesian
- **Middle Eastern**: Arabic, Persian, Hebrew, Turkish

### Other Languages
- **African**: Afrikaans
- **Minority/regional**: Basque, Galician, Armenian, Georgian

See [Stanza's complete model list](https://stanfordnlp.github.io/stanza/available_models.html) for detailed language coverage.

## Model Components

Depending on the language and package, each install may include:
- **Tokenization**: Sentence and token segmentation
- **POS tagging**: Universal POS tags and morphological features
- **Lemmatization**: Base form reduction
- **Dependency parsing**: Universal Dependencies syntax
- **NER**: Named entity recognition (subset of languages)
- **Sentiment**: Per-sentence sentiment (subset of languages)
- **Constituency**: Phrase-structure parse trees (`default`/`default_accurate` only)

Models are pretrained on Universal Dependencies v2.12 treebanks.

## Installation Process

### Admin Setup
1. **Install this data manager**: `data_manager_stanza_models`
2. **Install the Stanza tool**: `stanza_nlp`
3. **Navigate to Admin → Local Data**
4. **Select "Stanza Language Models"**

### Model Installation
1. **Choose a package**: `default_fast` (default), `default`, or `default_accurate`
2. **Choose languages**: Select checkboxes for the desired languages
3. **Run installation**: The data manager downloads and installs the models
4. **Verify installation**: Models appear in the Stanza tool's language dropdown, labelled with their package (e.g. `English — default_fast`)

### Post-Installation
- Models are immediately available to the Stanza NLP tool
- No restart required
- Models persist across Galaxy restarts

## Data Table Format

Models are registered in `stanza_models.loc` with five columns:
```
<value>    <name>    <lang>    <package>    <models_path>
```

- `value`: unique identifier, `<lang>-<package>` (also the on-disk subdirectory name)
- `name`: display name shown in the tool UI
- `lang`: ISO 639-1 language code
- `package`: `default_fast`, `default`, or `default_accurate`
- `models_path`: path to the `stanza_resources` directory containing the model

Example:
```
en-default_fast    English — default_fast    en    default_fast    /galaxy/tool-data/stanza_models/en-default_fast
en-default         English — default         en    default         /galaxy/tool-data/stanza_models/en-default
```

## Technical Details

### Download Source
- **Repository**: https://huggingface.co/stanfordnlp/
- **Model naming**: `stanza-{lang}` (e.g., `stanza-en`, `stanza-de`)
- **Fetched via**: `stanza.download()`

### Storage Structure
```
tool-data/
└── stanza_models/
    ├── en-default_fast/
    │   └── [English default_fast model files]
    ├── en-default/
    │   └── [English default model files, incl. constituency]
    └── stanza_models.loc
```

### Dependencies
- **Python 3.12**
- **stanza 1.12.0**

## Troubleshooting

### Common Issues
- **Network connectivity**: Ensure access to huggingface.co
- **Disk space**: Large language sets and the `default`/`default_accurate` packages require substantial storage
- **Permissions**: Galaxy must have write access to the managed data directory

### Model Verification
- Check `stanza_models.loc` for registered models
- Verify model files exist in the expected directories
- Test with the Stanza NLP tool after installation

## Citation

This data manager installs models created by the Stanford NLP Group. Please cite:

```
Qi, Peng, Yuhao Zhang, Yuhui Zhang, Jason Bolton, and Christopher D. Manning.
"Stanza: A Python Natural Language Processing Toolkit for Many Human Languages."
In Proceedings of the 58th Annual Meeting of the Association for Computational
Linguistics: System Demonstrations, 2020.
```
