#!/usr/bin/env python
# Copyright 2006 The Galaxy Project. All rights reserved.
"""
spaCy NLP Processing Script for Galaxy

Processes text files with spaCy and outputs results in various formats.
"""

import argparse
import json

import spacy


def process_text(nlp, text, output_format, include_components):
    """
    Process text with spaCy and format output.

    Args:
        nlp: spaCy language model
        text: Input text to process
        output_format: Output format (json, conll, conllu, text)
        include_components: List of components to include in output

    Returns:
        Formatted output string
    """
    doc = nlp(text)

    if output_format == "json":
        return format_json(doc, include_components)
    elif output_format == "conll":
        return format_conll(doc)
    elif output_format == "conllu":
        return format_conllu(doc)
    elif output_format == "text":
        return format_text(doc, include_components)
    else:
        return format_json(doc, include_components)


def format_json(doc, include_components):
    """Format document as JSON."""
    output = {
        "text": doc.text,
        "tokens": []
    }

    for token in doc:
        token_data = {
            "text": token.text,
            "start": token.idx,
            "end": token.idx + len(token.text),
            "is_alpha": token.is_alpha,
            "is_stop": token.is_stop,
        }

        if "pos" in include_components:
            token_data["pos"] = token.pos_
            token_data["tag"] = token.tag_

        if "lemma" in include_components:
            token_data["lemma"] = token.lemma_

        if "parse" in include_components:
            token_data["dep"] = token.dep_
            token_data["head"] = token.head.i

        output["tokens"].append(token_data)

    if "ner" in include_components:
        output["entities"] = [
            {
                "text": ent.text,
                "label": ent.label_,
                "start": ent.start_char,
                "end": ent.end_char,
                "start_token": ent.start,
                "end_token": ent.end
            }
            for ent in doc.ents
        ]

    output["sentences"] = [
        {
            "text": sent.text,
            "start": sent.start_char,
            "end": sent.end_char,
            "start_token": sent.start,
            "end_token": sent.end
        }
        for sent in doc.sents
    ]

    return json.dumps(output, indent=2, ensure_ascii=False)


def format_conll(doc):
    """Format document as CoNLL (tab-separated)."""
    lines = []
    for sent in doc.sents:
        for token in sent:
            # CoNLL format: ID, FORM, LEMMA, POS, NER, HEAD, DEPREL
            head_idx = token.head.i - sent.start + 1 if token.head.i != token.i else 0
            ner_tag = token.ent_type_ if token.ent_type_ else "O"

            line = f"{token.i - sent.start + 1}\t{token.text}\t{token.lemma_}\t{token.tag_}\t{ner_tag}\t{head_idx}\t{token.dep_}"
            lines.append(line)
        lines.append("")  # Empty line between sentences

    return "\n".join(lines)


def format_conllu(doc):
    """Format document as CoNLL-U (Universal Dependencies format)."""
    lines = []
    for sent in doc.sents:
        for token in sent:
            # CoNLL-U format: ID, FORM, LEMMA, UPOS, XPOS, FEATS, HEAD, DEPREL, DEPS, MISC
            head_idx = token.head.i - sent.start + 1 if token.head.i != token.i else 0

            # Get morphological features if available
            morph = str(token.morph) if token.morph else "_"

            line = f"{token.i - sent.start + 1}\t{token.text}\t{token.lemma_}\t{token.pos_}\t{token.tag_}\t{morph}\t{head_idx}\t{token.dep_}\t_\t_"
            lines.append(line)
        lines.append("")  # Empty line between sentences

    return "\n".join(lines)


def format_text(doc, include_components):
    """Format document as human-readable text."""
    lines = []

    # Document statistics
    num_tokens = len(doc)
    num_sents = len(list(doc.sents))
    lines.append(f"Document Statistics: {num_sents} sentences, {num_tokens} tokens\n")

    # Sentence-level output
    for i, sent in enumerate(doc.sents, 1):
        lines.append(f"\nSentence #{i} ({len(sent)} tokens):")
        lines.append(sent.text)
        lines.append("")

        if "pos" in include_components or "lemma" in include_components or "parse" in include_components:
            for token in sent:
                parts = [f"  {token.text}"]

                if "lemma" in include_components:
                    parts.append(f"lemma={token.lemma_}")

                if "pos" in include_components:
                    parts.append(f"pos={token.pos_}/{token.tag_}")

                if "parse" in include_components:
                    parts.append(f"dep={token.dep_}")
                    parts.append(f"head={token.head.text}")

                lines.append(" | ".join(parts))
            lines.append("")

    # Named entities
    if "ner" in include_components and doc.ents:
        lines.append("\nNamed Entities:")
        for ent in doc.ents:
            lines.append(f"  {ent.text} ({ent.label_})")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Process text with spaCy NLP")
    parser.add_argument("--input", required=True, help="Input text file")
    parser.add_argument("--output", required=True, help="Output file")
    parser.add_argument("--model-dir", required=True,
                        help="Path to a spaCy model directory (from the data manager)")
    parser.add_argument("--format", choices=["json", "conll", "conllu", "text"],
                        default="json", help="Output format")
    parser.add_argument("--annotators", required=True,
                        help="Comma-separated annotation types (tokenize, pos, ner, parse)")

    args = parser.parse_args()

    # Map annotator selection to components
    component_map = {
        "tokenize": ["sentences"],
        "pos": ["pos", "lemma"],
        "ner": ["pos", "lemma", "ner"],
        "parse": ["pos", "lemma", "parse"],
    }

    annotators = [a.strip() for a in args.annotators.split(",") if a.strip()]
    invalid = [a for a in annotators if a not in component_map]
    if invalid:
        parser.error(f"unknown annotation type(s): {', '.join(invalid)}")

    # Union the components of all selected annotators. Tokenization and sentence
    # segmentation are always included.
    selected = set(annotators) or {"tokenize"}
    include_components = list(
        {c for annotator in selected for c in component_map[annotator]}
    )

    # Load the data-manager-installed model directory
    nlp = spacy.load(args.model_dir)

    # Read input text
    with open(args.input, 'r', encoding='utf-8') as f:
        text = f.read()

    # Process text
    output = process_text(nlp, text, args.format, include_components)

    # Write output
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(output)

    print(f"Successfully processed {len(text)} characters")


if __name__ == "__main__":
    main()
