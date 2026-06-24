# Stanford CoreNLP Language Models Data Manager

Galaxy data manager for downloading and installing Stanford CoreNLP language model JARs.

## Description

This data manager allows Galaxy administrators to easily download and install Stanford CoreNLP language model JAR files from Maven Central for use with the Stanford CoreNLP annotation tool.

## Features

- **Multi-select support**: Install multiple language models at once using checkboxes
- **Common models option**: Optionally install common models JAR (required for coreference resolution)
- **Automatic download**: Downloads model JARs from Maven Central
- **Automatic registration**: Registers models in Galaxy's data table for immediate use
- **Error handling**: Continues with remaining models if one fails to download

## Supported Languages

Language models are available for:
- Arabic (ar)
- Chinese (zh)
- English (en)
- French (fr)
- German (de)
- Hungarian (hu)
- Italian (it)
- Spanish (es)

**Note**: Not all annotators are available for all languages. See the [Stanford CoreNLP documentation](https://stanfordnlp.github.io/CoreNLP/human-languages.html) for language-specific capabilities.

## Installation

1. Install this data manager via the Galaxy Tool Shed or manually
2. As a Galaxy admin, go to **Admin → Local Data**
3. Select "Stanford CoreNLP Language Models"
4. Check the boxes for the languages you want to install
5. Optionally check "Install common models" if you need coreference resolution (checked by default)
6. Click "Execute"

The models are large files (typically 100-500 MB each), so downloading multiple models may take several minutes depending on your connection speed.

### Common Models

The common models JAR (452 MB) contains shared dictionaries and resources needed for coreference resolution. If you plan to use the coreference annotator, you must install the common models. This option is checked by default for convenience.

## Requirements

- Python 3.9+
- Internet connection for downloading models from Maven Central

## Version

This data manager downloads models for CoreNLP version 4.5.10.

## Usage with Stanford CoreNLP Tool

After installing language models via this data manager, they will be available in the Stanford CoreNLP tool's "Language Model" dropdown. The tool uses a Docker container (`ksuderman/corenlp:4.5.10`) that includes the base CoreNLP library, and mounts the language-specific model JARs at runtime.

## More Information

- Stanford CoreNLP website: https://stanfordnlp.github.io/CoreNLP/
- Maven Central repository: https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/
- Galaxy Tool Shed: https://toolshed.g2.bx.psu.edu
