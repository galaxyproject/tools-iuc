#!/usr/bin/env python
# Copyright 2006 The Galaxy Project. All rights reserved.
"""
Stanza NLP Processing Script for Galaxy

Processes text files with Stanza and outputs results in various formats.

Author: Keith Suderman
License: MIT
"""

import argparse
import json
import os
import sys

try:
    import stanza
except ImportError:
    print("Error: Stanza is not installed. Please install stanza.", file=sys.stderr)
    sys.exit(1)


# Map annotator selections to Stanza processor strings
PROCESSOR_MAP = {
    "tokenize": "tokenize",
    "pos": "tokenize,mwt,pos,lemma",
    "ner": "tokenize,mwt,ner",
    "parse": "tokenize,mwt,pos,lemma,depparse",
    "sentiment": "tokenize,mwt,sentiment",
    "constituency": "tokenize,mwt,pos,constituency",
}


def process_text(doc, output_format, annotator):
    """Process a Stanza Document and format output."""
    if output_format == "json":
        return format_json(doc, annotator)
    elif output_format == "conll":
        return format_conll(doc)
    elif output_format == "conllu":
        return format_conllu(doc)
    elif output_format == "text":
        return format_text(doc, annotator)
    else:
        return format_json(doc, annotator)


def format_json(doc, annotator):
    """Format document as JSON."""
    output = {"text": doc.text, "sentences": []}

    for sent in doc.sentences:
        sent_data = {"text": sent.text, "tokens": []}

        for word in sent.words:
            token_data = {
                "text": word.text,
                "start_char": word.start_char,
                "end_char": word.end_char,
            }

            if annotator in ("pos", "parse", "constituency"):
                token_data["upos"] = word.upos
                token_data["xpos"] = word.xpos
                token_data["lemma"] = word.lemma
                if word.feats:
                    token_data["feats"] = word.feats

            if annotator == "parse":
                token_data["deprel"] = word.deprel
                token_data["head"] = word.head

            sent_data["tokens"].append(token_data)

        if annotator == "ner" and sent.ents:
            sent_data["entities"] = [
                {
                    "text": ent.text,
                    "type": ent.type,
                    "start_char": ent.start_char,
                    "end_char": ent.end_char,
                }
                for ent in sent.ents
            ]

        if annotator == "sentiment" and sent.sentiment is not None:
            sent_data["sentiment"] = sent.sentiment

        if annotator == "constituency" and sent.constituency is not None:
            sent_data["constituency"] = str(sent.constituency)

        output["sentences"].append(sent_data)

    return json.dumps(output, indent=2, ensure_ascii=False)


def format_conll(doc):
    """Format document as CoNLL (tab-separated)."""
    lines = []
    for sent in doc.sentences:
        for word in sent.words:
            ner_tag = "O"
            if hasattr(word, 'parent') and word.parent and hasattr(word.parent, 'ner'):
                ner_tag = word.parent.ner if word.parent.ner else "O"
            head = word.head if word.head is not None else 0
            deprel = word.deprel if word.deprel else "_"
            lemma = word.lemma if word.lemma else "_"
            xpos = word.xpos if word.xpos else "_"

            line = f"{word.id}\t{word.text}\t{lemma}\t{xpos}\t{ner_tag}\t{head}\t{deprel}"
            lines.append(line)
        lines.append("")
    return "\n".join(lines)


def format_conllu(doc):
    """Format document as CoNLL-U (Universal Dependencies format)."""
    lines = []
    for sent in doc.sentences:
        for word in sent.words:
            upos = word.upos if word.upos else "_"
            xpos = word.xpos if word.xpos else "_"
            lemma = word.lemma if word.lemma else "_"
            feats = word.feats if word.feats else "_"
            head = word.head if word.head is not None else 0
            deprel = word.deprel if word.deprel else "_"

            line = f"{word.id}\t{word.text}\t{lemma}\t{upos}\t{xpos}\t{feats}\t{head}\t{deprel}\t_\t_"
            lines.append(line)
        lines.append("")
    return "\n".join(lines)


def format_text(doc, annotator):
    """Format document as human-readable text."""
    lines = []

    num_tokens = sum(len(sent.words) for sent in doc.sentences)
    num_sents = len(doc.sentences)
    lines.append(f"Document Statistics: {num_sents} sentences, {num_tokens} tokens\n")

    for i, sent in enumerate(doc.sentences, 1):
        lines.append(f"\nSentence #{i} ({len(sent.words)} tokens):")
        lines.append(sent.text)
        lines.append("")

        if annotator in ("pos", "parse", "constituency"):
            for word in sent.words:
                parts = [f"  {word.text}"]
                parts.append(f"lemma={word.lemma}")
                parts.append(f"upos={word.upos}")
                if word.xpos:
                    parts.append(f"xpos={word.xpos}")
                if annotator == "parse" and word.deprel:
                    parts.append(f"deprel={word.deprel}")
                    parts.append(f"head={word.head}")
                lines.append(" | ".join(parts))
            lines.append("")

        if annotator == "ner" and sent.ents:
            lines.append("  Named Entities:")
            for ent in sent.ents:
                lines.append(f"    {ent.text} ({ent.type})")
            lines.append("")

        if annotator == "sentiment" and sent.sentiment is not None:
            labels = {0: "negative", 1: "neutral", 2: "positive"}
            lines.append(f"  Sentiment: {labels.get(sent.sentiment, sent.sentiment)}")
            lines.append("")

        if annotator == "constituency" and sent.constituency is not None:
            lines.append(f"  Constituency: {sent.constituency}")
            lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Process text with Stanza NLP")
    parser.add_argument("--input", required=True, help="Input text file")
    parser.add_argument("--output", required=True, help="Output file")
    parser.add_argument("--lang", required=True, help="Language code")
    parser.add_argument("--model-dir", required=True, help="Path to stanza_resources directory")
    parser.add_argument("--format", choices=["json", "conll", "conllu", "text"],
                        default="json", help="Output format")
    parser.add_argument("--annotators", required=True, help="Annotation type")
    parser.add_argument("--package", default="default_fast",
                        help="Stanza model package to load (e.g. default_fast, default, "
                             "default_accurate). Determines which annotators are available; "
                             "constituency parsing requires default or default_accurate.")
    parser.add_argument("--download", action="store_true",
                        help="Download the required models on the fly before processing. "
                             "Intended for testing; production models are installed via the data manager.")

    args = parser.parse_args()

    processors = PROCESSOR_MAP.get(args.annotators, "tokenize")

    # The processor that must be present for the selected annotation type. If the
    # chosen package does not ship this model, Stanza silently drops it, so we
    # check for it after loading and fail with a clear message instead.
    REQUIRED_PROCESSOR = {
        "tokenize": "tokenize",
        "pos": "pos",
        "ner": "ner",
        "parse": "depparse",
        "sentiment": "sentiment",
        "constituency": "constituency",
    }

    # In download mode, fetch only the models needed for the selected processors into
    # a writable directory in the job working area. This keeps tests self-contained
    # without bundling model files or relying on data-manager-installed models.
    model_dir = args.model_dir
    if args.download:
        model_dir = os.path.join(os.getcwd(), "stanza_resources")
        try:
            stanza.download(
                lang=args.lang,
                model_dir=model_dir,
                processors=processors,
                package=args.package,
                verbose=False,
            )
        except Exception as e:
            print(f"Error downloading Stanza model: {e}", file=sys.stderr)
            sys.exit(1)

    # Load Stanza pipeline using the selected package. default_fast (nocharlm) is
    # the low-memory default; default/default_accurate add heavier models such as
    # constituency parsing.
    try:
        nlp = stanza.Pipeline(
            lang=args.lang,
            dir=model_dir,
            processors=processors,
            package=args.package,
            download_method=None,
            use_gpu=False,
        )
    except Exception as e:
        print(f"Error loading Stanza pipeline: {e}", file=sys.stderr)
        sys.exit(1)

    # Stanza silently ignores a requested processor the package does not provide
    # (e.g. constituency is absent from default_fast). Fail clearly rather than
    # emitting empty annotations.
    required = REQUIRED_PROCESSOR.get(args.annotators)
    if required and required not in nlp.processors:
        print(
            f"Error: the '{args.annotators}' annotation requires the '{required}' model, "
            f"which is not available in the '{args.package}' package for language "
            f"'{args.lang}'. Install a package that includes it (default or "
            f"default_accurate) via the Stanza Language Models data manager.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Read input text
    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            text = f.read()
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Process text
    try:
        doc = nlp(text)
    except Exception as e:
        print(f"Error processing text: {e}", file=sys.stderr)
        sys.exit(1)

    # Format and write output
    try:
        output = process_text(doc, args.format, args.annotators)
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
    except Exception as e:
        print(f"Error writing output: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully processed {len(text)} characters")


if __name__ == "__main__":
    main()
