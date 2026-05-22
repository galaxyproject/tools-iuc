#!/usr/bin/env python
"""
VADER Sentiment Analysis for Galaxy

Performs sentence-level or document-level sentiment analysis using VADER.
Uses the conda vadersentiment package.
"""

import argparse
import csv
import json
import os
import re
import sys

from vaderSentiment.vaderSentiment import SentimentIntensityAnalyzer

def split_sentences(text):
    """Split text into sentences using a simple regex-based approach."""
    # Split on sentence-ending punctuation followed by whitespace
    raw = re.split(r'(?<=[.!?])\s+', text.strip())
    # Filter out empty strings and whitespace-only
    return [s.strip() for s in raw if s.strip()]


def classify(compound):
    """Classify compound score as positive/negative/neutral."""
    if compound >= 0.05:
        return "positive"
    elif compound <= -0.05:
        return "negative"
    else:
        return "neutral"


def analyze_sentences(text, analyzer):
    """Analyze each sentence individually."""
    sentences = split_sentences(text)
    results = []
    for sent in sentences:
        scores = analyzer.polarity_scores(sent)
        results.append({
            "text": sent,
            "compound": scores["compound"],
            "positive": scores["pos"],
            "negative": scores["neg"],
            "neutral": scores["neu"],
            "label": classify(scores["compound"]),
        })
    return results


def analyze_document(text, analyzer):
    """Analyze the entire document as one unit."""
    scores = analyzer.polarity_scores(text)
    return [{
        "text": text[:200] + ("..." if len(text) > 200 else ""),
        "compound": scores["compound"],
        "positive": scores["pos"],
        "negative": scores["neg"],
        "neutral": scores["neu"],
        "label": classify(scores["compound"]),
    }]


def write_tsv(results, output_path):
    """Write results as TSV."""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["sentence", "compound", "positive", "negative", "neutral", "label"])
        for r in results:
            writer.writerow([
                r["text"],
                f"{r['compound']:+.4f}",
                f"{r['positive']:.4f}",
                f"{r['negative']:.4f}",
                f"{r['neutral']:.4f}",
                r["label"],
            ])


def write_json(results, output_path):
    """Write results as JSON."""
    output = {
        "sentences": results,
        "summary": {
            "total": len(results),
            "positive": sum(1 for r in results if r["label"] == "positive"),
            "negative": sum(1 for r in results if r["label"] == "negative"),
            "neutral": sum(1 for r in results if r["label"] == "neutral"),
            "mean_compound": sum(r["compound"] for r in results) / len(results) if results else 0,
        }
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser(description="VADER sentiment analysis")
    parser.add_argument("--input", required=True, help="Input text file")
    parser.add_argument("--output", required=True, help="Output file")
    parser.add_argument("--format", choices=["tsv", "json"], default="tsv")
    parser.add_argument("--granularity", choices=["sentence", "document"], default="sentence")

    args = parser.parse_args()

    # Read input
    with open(args.input, 'r', encoding='utf-8') as f:
        text = f.read()

    if not text.strip():
        print("Warning: empty input file", file=sys.stderr)
        results = []
    else:
        # Initialize VADER with bundled lexicon
        lexicon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vader_lexicon.txt')
        analyzer = SentimentIntensityAnalyzer(lexicon_file=lexicon_path)

        if args.granularity == "sentence":
            results = analyze_sentences(text, analyzer)
        else:
            results = analyze_document(text, analyzer)

    # Write output
    if args.format == "tsv":
        write_tsv(results, args.output)
    else:
        write_json(results, args.output)

    # Summary
    if results:
        pos = sum(1 for r in results if r["label"] == "positive")
        neg = sum(1 for r in results if r["label"] == "negative")
        neu = sum(1 for r in results if r["label"] == "neutral")
        mean = sum(r["compound"] for r in results) / len(results)
        print(f"Analyzed {len(results)} unit(s): {pos} positive, {neg} negative, {neu} neutral")
        print(f"Mean compound score: {mean:+.4f}")


if __name__ == "__main__":
    main()
