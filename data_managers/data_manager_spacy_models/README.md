# spaCy Language Models Data Manager

Galaxy data manager for downloading and installing spaCy language models.

## Description

This data manager allows Galaxy administrators to easily download and install spaCy language models for use with the spaCy NLP tool. It supports 70+ models across 25+ languages.

## Features

- **Multi-select support**: Install multiple language models at once using checkboxes
- **Automatic installation**: Downloads models via `python -m spacy download`
- **Automatic registration**: Registers models in Galaxy's data table for immediate use
- **Error handling**: Continues with remaining models if one fails to download

## Supported Models

Models are available for 25+ languages including:
- English, Spanish, German, French
- Chinese, Japanese, Korean
- Portuguese, Italian, Dutch
- Greek, Polish, Norwegian
- And many more...

Each language typically has multiple model sizes:
- **sm** (small): Fast, less accurate, ~10-20MB
- **md** (medium): Balanced, includes word vectors, ~40-50MB
- **lg** (large): Most accurate, large word vectors, ~500MB-800MB
- **trf** (transformer): Transformer-based, highest accuracy, ~400MB-500MB

## Installation

1. Install this data manager via the Galaxy Tool Shed or manually
2. As a Galaxy admin, go to **Admin → Local Data**
3. Select "spaCy Language Models"
4. Check the boxes for the models you want to install
5. Click "Execute"

## Requirements

- spaCy 3.7.0 (automatically installed via conda)
- Internet connection for downloading models

## Version

This data manager is compatible with spaCy version 3.7.0.

## More Information

- spaCy website: https://spacy.io
- spaCy models documentation: https://spacy.io/usage/models
- Galaxy Tool Shed: https://toolshed.g2.bx.psu.edu
