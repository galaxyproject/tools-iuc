# Galaxy Wrapper for VADER Sentiment Analysis

This Galaxy tool performs sentiment analysis using VADER (Valence Aware Dictionary and sEntiment Reasoner), a lexicon and rule-based sentiment analysis tool specifically designed for social media and informal text.

## Features

- **Social media optimized**: Handles slang, emoticons, punctuation emphasis, and informal language
- **No dependencies**: Bundled lexicon with no external model downloads required
- **Fast processing**: Rule-based approach for high-speed sentiment analysis
- **Interpretable scores**: Clear numerical scores with rule-based explanations
- **Flexible granularity**: Per-sentence or whole-document analysis
- **Dual output formats**: TSV for analysis and JSON for programmatic use

## Requirements

- **Input**: Plain text files
- **No dependencies**: Pure Python implementation with bundled VADER lexicon
- **No setup required**: Ready to use immediately after installation

## When to Use VADER

VADER excels at analyzing:

| Text Type | Why VADER Works Well |
|---|---|
| **Political speeches** | Handles rhetorical language, emphasis, and persuasive tone |
| **News articles** | Captures positive/negative framing and editorial stance |
| **Social media posts** | Originally designed for Twitter, Facebook, and informal text |
| **Opinion pieces** | Effective on subjective, evaluative text |
| **Short informal text** | Handles slang, emoticons, and unconventional punctuation |

## VADER vs. Model-Based Approaches

Choose VADER over neural sentiment models (spaCy, Stanza) when:

- ✅ **Interpretability**: Need rule-based, explainable sentiment scores
- ✅ **Speed**: Processing large volumes of text quickly
- ✅ **No models**: Don't want to download large language models
- ✅ **Social media**: Analyzing informal, social media-style text
- ✅ **Lightweight**: Minimal computational requirements

Choose neural models for:
- 📚 **Formal text**: Academic papers, literature, formal documents
- 🎭 **Complex sentiment**: Subtle, contextual, or sarcastic expressions
- 🌍 **Multilingual**: Non-English text analysis

## Sentiment Scores

VADER produces four scores for each text unit:

### Compound Score (-1 to +1)
**Primary metric**: Overall sentiment intensity
- `≥ +0.05`: **Positive** sentiment
- `≤ -0.05`: **Negative** sentiment  
- `-0.05 to +0.05`: **Neutral** sentiment

### Component Scores (0 to 1)
- **Positive**: Proportion of positive sentiment
- **Negative**: Proportion of negative sentiment
- **Neutral**: Proportion of neutral sentiment

*Note: Positive + Negative + Neutral = 1.0*

## Analysis Granularity

### Per Sentence (Default)
- Splits text into sentences automatically
- Scores each sentence individually
- Produces table with one row per sentence
- Best for: Tracking sentiment changes throughout a document

### Whole Document
- Scores entire text as single unit
- Single sentiment score for the complete document  
- Best for: Overall document sentiment classification

## Output Formats

### Tabular (TSV)
Tab-separated file with columns:
- `sentence`: Text of the analyzed unit
- `compound`: Overall sentiment score (-1 to +1)
- `positive`: Positive proportion (0 to 1)
- `negative`: Negative proportion (0 to 1)
- `neutral`: Neutral proportion (0 to 1)
- `label`: Classification (positive/negative/neutral)

### JSON
Structured JSON format:
```json
[
  {
    "text": "This movie is great!",
    "compound": 0.6249,
    "positive": 0.779,
    "negative": 0.0,
    "neutral": 0.221,
    "label": "positive"
  }
]
```

## Example Use Cases

- **Political analysis**: Sentiment tracking in campaign speeches and debates
- **Media monitoring**: News sentiment analysis and bias detection
- **Social media research**: Public opinion analysis on Twitter/Facebook
- **Customer feedback**: Product review sentiment classification
- **Historical text analysis**: Sentiment trends in historical documents
- **Content moderation**: Identifying negative sentiment in user-generated content

## VADER's Strengths

- **Punctuation awareness**: "Good!!!" is more positive than "Good"
- **Capitalization**: "AMAZING" is stronger than "amazing"
- **Degree modifiers**: "very good" vs "slightly good"
- **Conjunction handling**: "but" and "however" shift sentiment
- **Emoticon support**: :) :( :D are properly interpreted
- **Slang recognition**: Modern informal expressions

## Installation

Install this tool from the Galaxy Toolshed: `vader_sentiment`

No additional setup required - the VADER lexicon is bundled with the tool.

## Citation

If you use this tool, please cite the original VADER paper:

```
Hutto, C.J. & Gilbert, Eric. (2014). VADER: A Parsimonious Rule-based Model for 
Sentiment Analysis of Social Media Text. Proceedings of the Eighth International 
AAAI Conference on Weblogs and Social Media.
```

## License

VADER is released under the MIT License, allowing unrestricted use in research and commercial applications.

## Version History

- **3.3.2+galaxy0**: Initial Galaxy release with bundled lexicon and dual output formats