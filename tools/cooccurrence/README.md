# Galaxy Wrapper for Co-occurrence Analysis

This Galaxy tool analyzes word co-occurrence relationships from NLP JSON output, enabling researchers to discover semantic and syntactic patterns in text corpora.

## Features

- **Multiple analysis methods**: Sentence-level, sliding window, and dependency-based co-occurrence detection
- **Flexible input**: Works with JSON output from spaCy, Stanza, or CoreNLP tools
- **Term representation options**: Lemma (recommended), surface form, or lowercased text
- **POS filtering**: Optional part-of-speech filtering (NOUN, PROPN, VERB, ADJ, ADV, NUM)
- **Stop word handling**: Built-in stop word removal plus support for custom stop word lists
- **Named entity analysis**: Option to restrict analysis to named entity co-occurrences only
- **Multiple outputs**: Co-occurrence pair list (TSV) and optional symmetric matrix
- **Configurable parameters**: Adjustable window sizes, frequency thresholds, and filtering options

## Requirements

- **Input**: JSON output from Galaxy NLP tools (spaCy, Stanza, or CoreNLP)
- **No dependencies**: Pure Python implementation with no external model downloads required

## Analysis Methods

| Method | Description | Use Case |
|---|---|---|
| **Sentence-level** | Terms co-occur if they appear in the same sentence | Document-level topic analysis (recommended starting point) |
| **Sliding window** | Terms co-occur within a fixed token window | Local semantic relationships and collocations |
| **Dependency-based** | Terms co-occur if connected by syntactic dependencies | Grammatical relationships (requires dependency parse) |

## Input Format

The tool expects JSON input with this structure from spaCy or Stanza:
```json
{
  "sentences": [
    {
      "tokens": [
        {
          "text": "word",
          "lemma": "lemmatized_form", 
          "pos": "POS_TAG",
          "is_alpha": true,
          "is_stop": false
        }
      ],
      "dependencies": [...] // for dependency-based analysis
    }
  ]
}
```

## Output Formats

### Co-occurrence Pairs (TSV)
Tab-separated file with columns:
- `term1`: First term in the pair
- `term2`: Second term in the pair  
- `count`: Number of co-occurrences

Results are sorted by count in descending order.

### Co-occurrence Matrix (TSV) 
Optional full term-by-term matrix where:
- Rows and columns represent vocabulary terms
- Cell values represent co-occurrence counts
- Can be large for extensive vocabularies

## Key Parameters

### Term Representation
- **Lemma** (recommended): Reduces inflected forms to base form ("supports" → "support")
- **Surface form**: Uses original text as-is
- **Lowercased**: Simple case normalization

### Filtering Options
- **POS tag restriction**: Focus on specific parts of speech (nouns, verbs, etc.)
- **Remove stop words**: Uses `is_stop` field from spaCy (for Stanza, use POS filtering)
- **Alphabetic only**: Exclude punctuation and numbers
- **Named entities only**: Restrict to named entity spans (PERSON, ORG, GPE, etc.)
- **Custom stop words**: Upload text file with one word per line
- **Minimum count**: Exclude pairs below frequency threshold

### Method-Specific Options
- **Window size**: For sliding window analysis (2-50 tokens, default: 5)

## Example Use Cases

- **Literary analysis**: Discover character relationships and thematic connections
- **Historical research**: Track concept associations across time periods  
- **Entity networks**: Build networks of people, organizations, and places
- **Corpus linguistics**: Identify collocation patterns and semantic fields
- **Digital humanities**: Analyze term associations in historical documents

## Example Workflow

1. Upload text → spaCy NLP (POS annotation, JSON output)
2. spaCy JSON → Co-occurrence Analysis (sentence-level, NOUN + PROPN, remove stops)
3. Pair list → downstream visualization or network analysis

## Installation

Install this tool from the Galaxy Toolshed: `cooccurrence_analysis`

No additional setup required - the tool is ready to use after installation.

## Citation

This tool implements foundational co-occurrence analysis methods. Please cite:

```
Manning, Christopher D. and Hinrich Schütze. 
Foundations of Statistical Natural Language Processing. 
MIT Press, 1999.
```

## Version History

- **1.0.0+galaxy0**: Initial release with sentence, window, and dependency-based analysis methods