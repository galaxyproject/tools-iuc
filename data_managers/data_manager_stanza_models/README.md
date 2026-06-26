# Galaxy Data Manager for Stanza Language Models

This Galaxy data manager downloads and installs Stanza language models for use with the Stanza NLP annotation tool, supporting 80+ languages with neural models trained on Universal Dependencies.

## Features

- **80+ languages**: Comprehensive language support for multilingual NLP
- **Direct HuggingFace download**: Downloads models directly from HuggingFace without requiring stanza installation
- **Multiple language installation**: Select and install multiple languages simultaneously
- **Progress reporting**: Shows download progress for each language model
- **Duplicate prevention**: Checks existing installations to avoid redundant downloads
- **Data table integration**: Automatically registers models in Galaxy's data table system

## How It Works

This data manager:
1. **Connects to HuggingFace**: Downloads default_fast model packages directly from Stanford's HuggingFace repository
2. **No dependencies**: Uses only Python's `urllib.request` - no stanza installation required
3. **Extracts models**: Unzips model packages to Galaxy's managed storage
4. **Registers models**: Updates the `stanza_models.loc` data table for tool access
5. **Version control**: Downloads models compatible with Stanza 1.11.1

## Supported Languages

The data manager supports **80+ languages** including:

### European Languages
- **Western**: English, German, French, Spanish, Italian, Portuguese, Dutch
- **Nordic**: Swedish, Danish, Norwegian (Bokmål/Nynorsk), Finnish
- **Slavic**: Russian, Ukrainian, Polish, Czech, Slovak, Croatian, Serbian, Bulgarian
- **Other**: Greek, Hungarian, Romanian, Estonian, Latvian, Lithuanian

### Asian Languages
- **East Asian**: Chinese (Simplified/Traditional), Japanese, Korean
- **South Asian**: Hindi, Tamil, Telugu, Marathi, Urdu
- **Southeast Asian**: Vietnamese, Thai, Indonesian
- **Middle Eastern**: Arabic, Persian, Hebrew, Turkish

### Other Languages
- **African**: Afrikaans
- **Minority**: Basque, Galician, Catalan, Armenian, Georgian

See [Stanza's complete model list](https://stanfordnlp.github.io/stanza/available_models.html) for detailed language coverage.

## Model Details

### Model Type
- **default_fast**: Memory-efficient models without character-level processing
- **Neural networks**: Pretrained on Universal Dependencies v2.12 treebanks
- **Multi-task**: Single package includes tokenization, POS, lemma, parsing, and NER models (where available)

### Model Sizes
- **Typical size**: 50-200MB per language
- **Variation**: Depends on language complexity and available training data
- **Storage**: Models persist in Galaxy's `tool-data/stanza_models/` directory

### Model Components
Each language package may include:
- **Tokenization**: Sentence and token segmentation
- **POS tagging**: Universal POS tags and morphological features
- **Lemmatization**: Base form reduction
- **Dependency parsing**: Universal Dependencies syntax
- **NER**: Named entity recognition (available for subset of languages)

## Installation Process

### Admin Setup
1. **Install this data manager**: `data_manager_stanza_models`
2. **Install the Stanza tool**: `stanza_nlp`
3. **Navigate to Admin → Local Data**
4. **Select "Stanza Language Models"**

### Model Installation
1. **Choose languages**: Select checkboxes for desired languages
2. **Run installation**: Data manager will download and extract models
3. **Monitor progress**: Download status shown for each language
4. **Verify installation**: Models appear in the Stanza tool's language dropdown

### Post-Installation
- Models are immediately available to the Stanza NLP tool
- No restart required
- Models persist across Galaxy restarts
- Multiple installations of the same language are prevented

## Data Table Format

Models are registered in `stanza_models.loc` with this format:
```
<lang_code>    <display_name>    <lang_code>    <models_path>
```

Example:
```
en    English    en    /galaxy/tool-data/stanza_models/en
de    German     de    /galaxy/tool-data/stanza_models/de
```

## Technical Details

### Download Source
- **Repository**: https://huggingface.co/stanfordnlp/
- **Model naming**: `stanza-{lang}` (e.g., `stanza-en`, `stanza-de`)
- **Version**: Models tagged with `v{resources_version}` from Stanford's resources.json

### Storage Structure
```
tool-data/
└── stanza_models/
    ├── en/
    │   └── [English model files]
    ├── de/
    │   └── [German model files]
    └── stanza_models.loc
```

### Dependencies
- **Python 3.12**: Standard library only
- **No stanza package**: Downloads directly from HuggingFace
- **urllib.request**: For HTTP downloads
- **zipfile**: For model extraction

## Troubleshooting

### Common Issues
- **Network connectivity**: Ensure access to huggingface.co
- **Disk space**: Large language sets require substantial storage
- **Permissions**: Galaxy must have write access to tool-data directory

### Model Verification
- Check `stanza_models.loc` for registered models
- Verify model files exist in expected directories
- Test with Stanza NLP tool after installation

## Citation

This data manager installs models created by the Stanford NLP Group. Please cite:

```
Qi, Peng, Yuhao Zhang, Yuhui Zhang, Jason Bolton, and Christopher D. Manning. 
"Stanza: A Python Natural Language Processing Toolkit for Many Human Languages." 
In Proceedings of the 58th Annual Meeting of the Association for Computational 
Linguistics: System Demonstrations, 2020.
```

## Version History

- **1.11.1.3**: Enhanced duplicate prevention and error handling
- **1.11.1.2**: Improved download progress reporting
- **1.11.1.1**: Direct HuggingFace download implementation
- **1.11.1.0**: Initial release