# Galaxy Wrapper for Stanford CoreNLP (Multi-language)

This Galaxy tool wrapper provides access to Stanford CoreNLP annotation pipelines with support for multiple languages. Users can select from various annotators and output formats.

See the [Stanford CoreNLP website](https://stanfordnlp.github.io/CoreNLP/) for more information about the underlying library.

## Features

- **Multi-language support**: Process text in Arabic, Chinese, English, French, German, Hungarian, Italian, and Spanish
- **Multiple annotators**: Tokenization, POS tagging, NER, dependency parsing, coreference, and sentiment analysis
- **Multiple output formats**: JSON, CoNLL, CoNLL-U, text, and XML
- **Dockerized execution**: Uses Docker container for consistent environment
- **Data manager integration**: Language models downloaded separately for flexibility

## Requirements

- **Data Manager**: Language models must be installed via the Stanford CoreNLP Language Models data manager
- **Common Models**: Required for coreference resolution (installed via data manager)
- **Docker**: The tool uses the `ksuderman/corenlp:4.5.10` Docker image

## Annotators

| Type | Description | Annotators Used |
|---|---|---|
| **Segmentation** | Sentences and tokens only | tokenize, ssplit |
| **Part of speech** | POS tags and lemmas | tokenize, ssplit, pos, lemma |
| **Named entity recognition** | Named entities (PERSON, ORG, etc.) | tokenize, ssplit, pos, lemma, ner |
| **Dependency parse** | Syntactic dependencies | tokenize, ssplit, pos, lemma, ner, parse |
| **Coreference** | Entity coreferences (requires common models) | tokenize, ssplit, pos, lemma, ner, parse, coref |
| **Sentiment analysis** | Sentiment scores | tokenize, ssplit, pos, lemma, parse, sentiment |

## Output Formats

- **JSON**: Full annotations with hierarchical structure (recommended for most use cases)
- **CoNLL**: Tab-separated format for NER and basic annotations
- **CoNLL-U**: Universal Dependencies format
- **Text**: Human-readable plain text output
- **XML**: Structured XML format

**Note**: Not all annotations are representable in all formats. JSON provides the most complete representation.

## Installation

1. Install the data manager: `data_manager_corenlp_models`
2. Install this tool: `stanford_corenlp_multilang`
3. Use the data manager to download language models:
   - Go to **Admin → Local Data**
   - Select "Stanford CoreNLP Language Models"
   - Choose language(s) to install
   - Check "Install common models" if you need coreference resolution
4. The tool will automatically appear in the "Natural Language Processing" section

## Version

This wrapper is for Stanford CoreNLP version 4.5.10.


## LICENSE

This tool wrapper is provided under the MIT License.