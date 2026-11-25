#!/usr/bin/env python
"""Build the tiny spaCy model used by the spaCy NLP tool tests.

The wrapper needs a loadable spaCy model directory in the ``spacy_models`` data
table. Real models (en_core_web_sm and up) are several MB, so the tests ship a
purpose-built model instead: a few kB of weights trained on the two sample
sentences below, plus a deterministic entity ruler and a lookup lemmatizer.

The model is not a copy of any published spaCy pipeline, so it carries no
third-party model licence, and its annotations are stable across spaCy patch
releases because the pipeline contains no statistical NER.

Regenerate with (spaCy 3.8, no extra packages needed):

    pip install "spacy>=3.8,<3.9"
    python test-data/build_test_model.py test-data/spacy_test_model

The output directory is what the ``spacy_models`` data table points at
(see ``test-data/spacy_models.loc``).
"""

import copy
import sys

import spacy
from spacy.lookups import Lookups
from spacy.training import Example

# Two short, fully annotated sentences. The tagger memorises their vocabulary
# (John, Smith, went, ...) so the tests get stable POS tags and lemmas; named
# entities come from the entity ruler below.
TRAIN = [
    {
        "text": "John Smith went to Walmart .",
        "tags": ["PROPN", "PROPN", "VERB", "ADP", "PROPN", "PUNCT"],
        "heads": [1, 2, 2, 2, 3, 2],
        "deps": ["compound", "nsubj", "ROOT", "prep", "pobj", "punct"],
    },
    {
        "text": "Mary Jones visited Berlin in May 2020 .",
        "tags": ["PROPN", "PROPN", "VERB", "PROPN", "ADP", "PROPN", "NUM", "PUNCT"],
        "heads": [1, 2, 2, 2, 2, 4, 5, 2],
        "deps": ["compound", "nsubj", "ROOT", "dobj", "prep", "pobj", "nummod", "punct"],
    },
]

PATTERNS = [
    {"label": "PERSON", "pattern": [{"LOWER": "john"}, {"LOWER": "smith"}]},
    {"label": "PERSON", "pattern": [{"LOWER": "mary"}, {"LOWER": "jones"}]},
    {"label": "ORG", "pattern": [{"LOWER": "walmart"}]},
    {"label": "ORG", "pattern": [{"LOWER": "ibm"}]},
    {"label": "GPE", "pattern": [{"LOWER": "berlin"}]},
    {"label": "DATE", "pattern": [
        {"LOWER": {"IN": ["january", "february", "march", "april", "may", "june", "july",
                          "august", "september", "october", "november", "december"]}},
        {"IS_DIGIT": True}]},
    {"label": "DATE", "pattern": [{"IS_DIGIT": True}, {"IS_PUNCT": True}, {"IS_DIGIT": True}]},
]

# A lemmatizer backed by the lookup table below: spacy-lookups-data is not a
# dependency of the container, and the table is saved inside the model.
LEMMAS = {
    "went": "go",
    "visited": "visit",
    "buy": "buy",
    "is": "be",
    "are": "be",
    "the": "the",
}


def build(nlp, n_iter=200):
    examples = []
    for rec in TRAIN:
        doc = nlp.make_doc(rec["text"])
        examples.append(Example.from_dict(doc, {
            "tags": rec["tags"],
            "pos": rec["tags"],
            "heads": rec["heads"],
            "deps": rec["deps"],
        }))
    optimizer = nlp.initialize(get_examples=lambda: examples)
    for _ in range(n_iter):
        for example in examples:
            nlp.update([example], sgd=optimizer, drop=0.0)
    return nlp


def main(out_dir):
    nlp = spacy.blank("en")
    # Hash embeddings keep the committed model well below 1 MB.
    tok2vec = {
        "@architectures": "spacy.HashEmbedCNN.v2",
        "width": 8,
        "depth": 2,
        "embed_size": 500,
        "window_size": 1,
        "maxout_pieces": 2,
        "subword_features": True,
        "pretrained_vectors": None,
    }
    nlp.add_pipe("tagger", config={"model": {"@architectures": "spacy.Tagger.v2",
                                             "tok2vec": copy.deepcopy(tok2vec)}})
    nlp.add_pipe("morphologizer", config={"model": {"@architectures": "spacy.Tagger.v2",
                                                    "tok2vec": copy.deepcopy(tok2vec)}})
    nlp.add_pipe("parser", config={
        # a handful of examples is all the parser gets
        "min_action_freq": 1,
        "model": {
            "@architectures": "spacy.TransitionBasedParser.v2",
            "state_type": "parser",
            "extra_state_tokens": False,
            "hidden_width": 8,
            "maxout_pieces": 2,
            "use_upper": True,
            "tok2vec": copy.deepcopy(tok2vec),
        }})

    build(nlp)

    nlp.add_pipe("entity_ruler").add_patterns(PATTERNS)

    lookups = Lookups()
    lookups.add_table("lemma_lookup", LEMMAS)
    nlp.add_pipe("lemmatizer", config={"mode": "lookup"}).initialize(nlp=nlp, lookups=lookups)

    nlp.meta.update({
        "name": "core_test_web_sm",
        "version": "1.0.0",
        "lang": "en",
        "description": "Tiny spaCy model built for the Galaxy spaCy wrapper tests",
        "author": "Galaxy IUC",
        "license": "MIT",
    })
    nlp.to_disk(out_dir)
    print(f"wrote {out_dir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "test-data/spacy_test_model")
