# Galaxy Wrapper for Stanford Stanza NLP

This Galaxy tool provides access to Stanza, Stanford's neural natural language processing toolkit, supporting 80+ human languages with state-of-the-art accuracy for multilingual text analysis.

## Features

- **80+ languages**: Comprehensive multilingual support for diverse text corpora
- **Neural models**: State-of-the-art accuracy with pretrained neural networks
- **Multiple annotators**: Tokenization, POS tagging, NER, parsing, sentiment, and constituency parsing
- **Universal Dependencies**: Standardized annotations following Universal Dependencies v2.12
- **Multiple output formats**: JSON, CoNLL, CoNLL-U, and human-readable text
- **Conda integration**: CPU-optimized PyTorch models via conda packages
- **Data manager integration**: Language models downloaded and managed separately

## Requirements

- **Data Manager**: Language models must be installed via the Stanza Language Models data manager
- **Python**: Python 3.12 with Stanza 1.12.0 and transformers 4.47.1 packages
- **Memory**: Uses default_fast models for efficient memory usage

## Annotation Types

| Annotator | Description | Output |
|---|---|---|
| **Tokenization** | Sentence segmentation and tokenization | Tokens with character offsets |
| **Part of speech** | POS tags, lemmas, and morphological features | Universal POS (UPOS), treebank POS (XPOS), lemmas |
| **Named entity recognition** | Person, organization, location, date entities | Entity spans with types (PERSON, ORG, GPE, etc.) |
| **Dependency parsing** | Syntactic dependencies following Universal Dependencies | Head-child relationships with dependency labels |
| **Sentiment analysis** | Per-sentence sentiment scoring | Sentiment scores (0=negative, 1=neutral, 2=positive) |
| **Constituency parsing** | Phrase structure parse trees | Hierarchical syntactic structure |

## Language Coverage

Stanza supports **80+ languages** including:

### Major Languages
- **European**: English, Spanish, German, French, Italian, Portuguese, Dutch, Swedish, Danish, Norwegian, Greek, Polish, Russian, Ukrainian
- **Asian**: Chinese, Japanese, Korean, Arabic, Hindi, Turkish  
- **Others**: And many more languages with Universal Dependencies treebanks

### NER Support
Named entity recognition is available for a subset of languages including:
- English, Chinese, Spanish, German, French, Dutch, Russian, Ukrainian

See [Stanza's model documentation](https://stanfordnlp.github.io/stanza/available_models.html) for the complete supported language list.

## Input Format

- **Text files**: Plain text input in any supported language
- **Encoding**: UTF-8 text encoding

## Output Formats

### JSON (Recommended)
Comprehensive structured output with all annotations:
```json
{
  "sentences": [
    {
      "tokens": [
        {
          "id": 1,
          "text": "John",
          "lemma": "John", 
          "upos": "PROPN",
          "head": 2,
          "deprel": "nsubj"
        }
      ],
      "entities": [
        {
          "text": "John Smith",
          "type": "PERSON",
          "start_char": 0,
          "end_char": 10
        }
      ]
    }
  ]
}
```

### CoNLL-U
Universal Dependencies format with morphological features:
```
1	John	John	PROPN	_	_	2	nsubj	_	_
2	works	work	VERB	_	_	0	root	_	_
```

### CoNLL
Tab-separated format suitable for dependency parsing analysis.

### Text
Human-readable output with statistics and formatted annotations.

## Model Architecture

- **Neural networks**: Pretrained neural models for each language and task
- **Universal Dependencies**: Consistent annotation standards across languages
- **Default-fast models**: Memory-efficient nocharlm models optimized for containers
- **CPU-optimized**: PyTorch models configured for CPU-only execution

## Example Use Cases

- **Multilingual corpus analysis**: Process text in 80+ languages with consistent annotations
- **Cross-lingual studies**: Compare linguistic phenomena across different languages
- **Historical linguistics**: Analyze texts in various languages and time periods
- **Digital humanities**: Multi-language support for international document collections
- **Dependency syntax**: Universal Dependencies parsing for computational linguistics

## Installation

1. Install the data manager: `data_manager_stanza_models`
2. Install this tool: `stanza_nlp`
3. Use the data manager to download language models:
   - Go to **Admin → Local Data**
   - Select "Stanza Language Models"
   - Choose language(s) to install
   - Models download directly from HuggingFace

## Performance Notes

- **Memory efficient**: Uses default_fast models without character-level modeling
- **CPU-optimized**: PyTorch configured for CPU-only execution
- **Conda environment**: Runs with conda dependency resolution for consistent package versions
- **Model caching**: Downloaded models persist across runs

## Citation

If you use this tool, please cite:

```
Qi, Peng, Yuhao Zhang, Yuhui Zhang, Jason Bolton, and Christopher D. Manning. 
"Stanza: A Python Natural Language Processing Toolkit for Many Human Languages." 
In Proceedings of the 58th Annual Meeting of the Association for Computational 
Linguistics: System Demonstrations, 2020.
```

## Testing

### CI Testing
The tool tests exercise real tokenization without bundling any model files:
- A hidden `download_models` parameter (enabled only in the tests) downloads the
  small `tokenize` and `mwt` models (~1 MB total) on the fly
- Only the models needed for the selected processors are fetched, so tests stay fast
- No large model files are committed to the repository, avoiding CI timeouts

### Production Usage
In production the `download_models` parameter stays off and language models are
installed via the data manager:
1. Use Admin → Local Data → "Stanza Language Models"
2. Download required language models from HuggingFace
3. Models will be cached for subsequent runs

## Version History

- **1.12.0+galaxy0**: Current release with updated Stanza version and lightweight testing
- **1.11.1+galaxy4**: Previous release with enhanced output formatting and CPU optimization
- **1.11.1+galaxy3**: Earlier stable release
- **1.11.1+galaxy2**: Earlier release
- **1.11.1+galaxy1**: Beta release