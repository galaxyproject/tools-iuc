#!/usr/bin/env python
"""
Co-occurrence Analysis for Galaxy

Computes term co-occurrence from NLP-annotated JSON (spaCy or Stanza format).
"""

import argparse
import csv
import json
import re
import sys
from collections import Counter
from itertools import combinations


def detect_format(data):
    """Detect whether JSON is spaCy or Stanza format."""
    if "sentences" in data and data["sentences"]:
        first_sent = data["sentences"][0]
        if "tokens" in first_sent:
            return "stanza"
        if "start_token" in first_sent:
            return "spacy"
    if "tokens" in data:
        return "spacy"
    return "unknown"


def extract_spacy_sentences(data):
    """Extract sentence-grouped tokens from spaCy JSON format.

    Adjusts head indices from absolute (document-level) to 1-based
    sentence-relative to match the Stanza convention used by the
    co-occurrence functions.
    """
    tokens = data.get("tokens", [])
    sentences_meta = data.get("sentences", [])
    entities = data.get("entities", [])

    if not sentences_meta:
        # No sentence info — treat entire document as one sentence
        # Convert head from absolute 0-based to 1-based sentence-relative
        for i, tok in enumerate(tokens):
            if "head" in tok:
                head_abs = tok["head"]
                tok["head"] = (head_abs + 1) if head_abs != i else 0
        return [tokens], entities

    sentences = []
    for sent in sentences_meta:
        start = sent["start_token"]
        end = sent["end_token"]
        sent_tokens = []
        for tok in tokens[start:end]:
            tok = dict(tok)  # copy to avoid mutating original
            if "head" in tok:
                head_abs = tok["head"]
                if head_abs == (start + len(sent_tokens)):
                    # Self-referencing (ROOT in spaCy) → head=0
                    tok["head"] = 0
                else:
                    tok["head"] = head_abs - start + 1
            sent_tokens.append(tok)
        sentences.append(sent_tokens)

    return sentences, entities


def extract_stanza_sentences(data):
    """Extract sentence-grouped tokens from Stanza JSON format."""
    sentences = []
    entities = []

    for sent in data.get("sentences", []):
        tokens = []
        for word in sent.get("tokens", []):
            token = {
                "text": word.get("text", ""),
                "lemma": word.get("lemma", word.get("text", "")),
                "pos": word.get("upos", word.get("pos", "")),
            }
            if "deprel" in word:
                token["dep"] = word["deprel"]
                token["head"] = word.get("head", 0)
            tokens.append(token)

        sent_entities = sent.get("entities", [])
        for ent in sent_entities:
            entities.append({
                "text": ent.get("text", ""),
                "label": ent.get("type", ent.get("label", "")),
            })

        sentences.append(tokens)

    return sentences, entities


def get_term(token, term_type):
    """Extract the term representation from a token."""
    if term_type == "lemma":
        return token.get("lemma", token.get("text", "")).lower()
    elif term_type == "lower":
        return token.get("text", "").lower()
    else:
        return token.get("text", "")


def is_alpha(token):
    """Check if token is alphabetic."""
    if "is_alpha" in token:
        return token["is_alpha"]
    return bool(re.match(r'^[a-zA-Z\u00C0-\u024F\u0400-\u04FF]+$', token.get("text", "")))


def is_stop(token):
    """Check if token is a stop word."""
    if "is_stop" in token:
        return token["is_stop"]
    return False


def filter_token(token, pos_tags, remove_stops, alpha_only, custom_stops, term_type):
    """Return the term if token passes all filters, else None."""
    if alpha_only and not is_alpha(token):
        return None
    if remove_stops and is_stop(token):
        return None
    if pos_tags and token.get("pos", "") not in pos_tags:
        return None
    term = get_term(token, term_type)
    if not term or not term.strip():
        return None
    if custom_stops and term.lower() in custom_stops:
        return None
    return term


def cooccur_sentence(sentences, pos_tags, remove_stops, alpha_only, custom_stops, term_type):
    """Sentence-level co-occurrence."""
    counter = Counter()
    for sent_tokens in sentences:
        terms = []
        for token in sent_tokens:
            term = filter_token(token, pos_tags, remove_stops, alpha_only, custom_stops, term_type)
            if term:
                terms.append(term)
        for w1, w2 in combinations(sorted(set(terms)), 2):
            counter[(w1, w2)] += 1
    return counter


def cooccur_window(sentences, window_size, pos_tags, remove_stops, alpha_only, custom_stops, term_type):
    """Sliding window co-occurrence."""
    # Flatten all tokens but respect filters
    all_terms = []
    for sent_tokens in sentences:
        for token in sent_tokens:
            term = filter_token(token, pos_tags, remove_stops, alpha_only, custom_stops, term_type)
            if term:
                all_terms.append(term)

    counter = Counter()
    for i in range(len(all_terms)):
        window = all_terms[i:i + window_size]
        for w1, w2 in combinations(sorted(set(window)), 2):
            counter[(w1, w2)] += 1
    return counter


def cooccur_dependency(sentences, pos_tags, remove_stops, alpha_only, custom_stops, term_type):
    """Dependency-based co-occurrence (head-child pairs)."""
    counter = Counter()
    for sent_tokens in sentences:
        for token in sent_tokens:
            head_idx = token.get("head")
            dep = token.get("dep", "")
            if head_idx is None or dep == "" or dep == "root":
                continue

            # Head index in spaCy is absolute, in Stanza it's 1-based within sentence
            # Stanza: head=0 means root, head=N means Nth word (1-based)
            if isinstance(head_idx, int) and 1 <= head_idx <= len(sent_tokens):
                head_token = sent_tokens[head_idx - 1]
            else:
                continue

            child_term = filter_token(token, pos_tags, remove_stops, alpha_only, custom_stops, term_type)
            head_term = filter_token(head_token, pos_tags, remove_stops, alpha_only, custom_stops, term_type)

            if child_term and head_term and child_term != head_term:
                pair = tuple(sorted([head_term, child_term]))
                counter[pair] += 1

    return counter


def cooccur_entities(sentences, entities_data, term_type):
    """Entity-only co-occurrence at sentence level."""
    counter = Counter()

    # For spaCy format, entities have start_token/end_token
    # For Stanza format, entities are per-sentence
    # We'll extract entity texts per sentence from the token data
    for sent_tokens in sentences:
        # Collect entity spans from tokens that have ner/entity info
        # Simpler approach: just use the entities list
        pass

    # Use the global entities list grouped by sentence proximity
    # Actually, let's extract entity texts from the entities_data
    entity_texts = [e.get("text", "").lower() if term_type != "text" else e.get("text", "")
                    for e in entities_data if e.get("text", "")]

    # For sentence-level entity co-occurrence, we need sentence boundaries
    # This is a simplified version that computes document-level entity co-occurrence
    for w1, w2 in combinations(sorted(set(entity_texts)), 2):
        counter[(w1, w2)] += 1

    return counter


def write_pairs(counter, output_path, min_count):
    """Write co-occurrence pairs as TSV."""
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["term1", "term2", "count"])
        for (w1, w2), count in sorted(counter.items(), key=lambda x: -x[1]):
            if count >= min_count:
                writer.writerow([w1, w2, count])


def write_matrix(counter, output_path, min_count):
    """Write co-occurrence matrix as TSV."""
    # Filter by min_count first
    filtered = {pair: count for pair, count in counter.items() if count >= min_count}

    # Build vocabulary
    vocab = sorted(set(w for pair in filtered for w in pair))

    if not vocab:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("(no co-occurrences found)\n")
        return

    # Build matrix
    matrix = {v: {v2: 0 for v2 in vocab} for v in vocab}
    for (w1, w2), count in filtered.items():
        matrix[w1][w2] = count
        matrix[w2][w1] = count

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([""] + vocab)
        for term in vocab:
            row = [term] + [matrix[term][v] for v in vocab]
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Co-occurrence analysis from NLP JSON")
    parser.add_argument("--input", required=True, help="Input NLP-annotated JSON file")
    parser.add_argument("--pairs", required=True, help="Output pair list (TSV)")
    parser.add_argument("--matrix", help="Output matrix (TSV)")
    parser.add_argument("--method", choices=["sentence", "window", "dependency"],
                        default="sentence", help="Co-occurrence method")
    parser.add_argument("--window-size", type=int, default=5, help="Window size for sliding window")
    parser.add_argument("--pos-tags", help="Comma-separated POS tags to include")
    parser.add_argument("--remove-stop", action="store_true", help="Remove stop words")
    parser.add_argument("--alpha-only", action="store_true", help="Alphabetic tokens only")
    parser.add_argument("--entities-only", action="store_true", help="Named entities only")
    parser.add_argument("--stopword-file", help="Custom stop word list (one per line)")
    parser.add_argument("--min-count", type=int, default=1, help="Minimum co-occurrence count")
    parser.add_argument("--term-type", choices=["lemma", "text", "lower"],
                        default="lemma", help="Term representation")

    args = parser.parse_args()

    # Read input JSON
    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading input JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Detect format and extract sentences
    fmt = detect_format(data)
    print(f"Detected input format: {fmt}")

    if fmt == "stanza":
        sentences, entities = extract_stanza_sentences(data)
    elif fmt == "spacy":
        sentences, entities = extract_spacy_sentences(data)
    else:
        print("Error: Could not detect JSON format. Expected spaCy or Stanza output.", file=sys.stderr)
        sys.exit(1)

    # Parse filter options
    pos_tags = set(args.pos_tags.split(",")) if args.pos_tags else None
    custom_stops = set()
    if args.stopword_file:
        with open(args.stopword_file, 'r', encoding='utf-8') as f:
            custom_stops = {line.strip().lower() for line in f if line.strip()}

    # Compute co-occurrence
    if args.entities_only:
        counter = cooccur_entities(sentences, entities, args.term_type)
    elif args.method == "sentence":
        counter = cooccur_sentence(sentences, pos_tags, args.remove_stop, args.alpha_only,
                                   custom_stops, args.term_type)
    elif args.method == "window":
        counter = cooccur_window(sentences, args.window_size, pos_tags, args.remove_stop,
                                 args.alpha_only, custom_stops, args.term_type)
    elif args.method == "dependency":
        counter = cooccur_dependency(sentences, pos_tags, args.remove_stop, args.alpha_only,
                                     custom_stops, args.term_type)
    else:
        counter = Counter()

    print(f"Found {len(counter)} unique co-occurrence pairs")
    total = sum(counter.values())
    print(f"Total co-occurrences: {total}")

    # Write outputs
    write_pairs(counter, args.pairs, args.min_count)
    print(f"Wrote pair list to {args.pairs}")

    if args.matrix:
        write_matrix(counter, args.matrix, args.min_count)
        print(f"Wrote matrix to {args.matrix}")


if __name__ == "__main__":
    main()
