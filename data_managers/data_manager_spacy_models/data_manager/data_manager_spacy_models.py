#!/usr/bin/env python
# Copyright 2006 The Galaxy Project. All rights reserved.
"""
Data Manager for spaCy Language Models

This script downloads and registers spaCy language models for use with Galaxy.
"""

import argparse
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

# GitHub releases base URL for spaCy models
SPACY_MODELS_BASE_URL = "https://github.com/explosion/spacy-models/releases/download"
# spaCy compatibility endpoint
SPACY_COMPAT_URL = "https://raw.githubusercontent.com/explosion/spacy-models/master/compatibility.json"
# spaCy version this data manager targets
SPACY_VERSION = "3.8"


# spaCy model information
# Format: model_name -> (display_name, language, size)
SPACY_MODELS = {
    # English
    "en_core_web_sm": ("English (small)", "en", "sm"),
    "en_core_web_md": ("English (medium)", "en", "md"),
    "en_core_web_lg": ("English (large)", "en", "lg"),
    "en_core_web_trf": ("English (transformer)", "en", "trf"),

    # Spanish
    "es_core_news_sm": ("Spanish (small)", "es", "sm"),
    "es_core_news_md": ("Spanish (medium)", "es", "md"),
    "es_core_news_lg": ("Spanish (large)", "es", "lg"),

    # German
    "de_core_news_sm": ("German (small)", "de", "sm"),
    "de_core_news_md": ("German (medium)", "de", "md"),
    "de_core_news_lg": ("German (large)", "de", "lg"),

    # French
    "fr_core_news_sm": ("French (small)", "fr", "sm"),
    "fr_core_news_md": ("French (medium)", "fr", "md"),
    "fr_core_news_lg": ("French (large)", "fr", "lg"),

    # Chinese
    "zh_core_web_sm": ("Chinese (small)", "zh", "sm"),
    "zh_core_web_md": ("Chinese (medium)", "zh", "md"),
    "zh_core_web_lg": ("Chinese (large)", "zh", "lg"),

    # Japanese
    "ja_core_news_sm": ("Japanese (small)", "ja", "sm"),
    "ja_core_news_md": ("Japanese (medium)", "ja", "md"),
    "ja_core_news_lg": ("Japanese (large)", "ja", "lg"),

    # Portuguese
    "pt_core_news_sm": ("Portuguese (small)", "pt", "sm"),
    "pt_core_news_md": ("Portuguese (medium)", "pt", "md"),
    "pt_core_news_lg": ("Portuguese (large)", "pt", "lg"),

    # Italian
    "it_core_news_sm": ("Italian (small)", "it", "sm"),
    "it_core_news_md": ("Italian (medium)", "it", "md"),
    "it_core_news_lg": ("Italian (large)", "it", "lg"),

    # Dutch
    "nl_core_news_sm": ("Dutch (small)", "nl", "sm"),
    "nl_core_news_md": ("Dutch (medium)", "nl", "md"),
    "nl_core_news_lg": ("Dutch (large)", "nl", "lg"),

    # Greek
    "el_core_news_sm": ("Greek (small)", "el", "sm"),
    "el_core_news_md": ("Greek (medium)", "el", "md"),
    "el_core_news_lg": ("Greek (large)", "el", "lg"),

    # Polish
    "pl_core_news_sm": ("Polish (small)", "pl", "sm"),
    "pl_core_news_md": ("Polish (medium)", "pl", "md"),
    "pl_core_news_lg": ("Polish (large)", "pl", "lg"),

    # Norwegian
    "nb_core_news_sm": ("Norwegian Bokmål (small)", "nb", "sm"),
    "nb_core_news_md": ("Norwegian Bokmål (medium)", "nb", "md"),
    "nb_core_news_lg": ("Norwegian Bokmål (large)", "nb", "lg"),

    # Lithuanian
    "lt_core_news_sm": ("Lithuanian (small)", "lt", "sm"),
    "lt_core_news_md": ("Lithuanian (medium)", "lt", "md"),
    "lt_core_news_lg": ("Lithuanian (large)", "lt", "lg"),

    # Danish
    "da_core_news_sm": ("Danish (small)", "da", "sm"),
    "da_core_news_md": ("Danish (medium)", "da", "md"),
    "da_core_news_lg": ("Danish (large)", "da", "lg"),

    # Swedish
    "sv_core_news_sm": ("Swedish (small)", "sv", "sm"),
    "sv_core_news_md": ("Swedish (medium)", "sv", "md"),
    "sv_core_news_lg": ("Swedish (large)", "sv", "lg"),

    # Romanian
    "ro_core_news_sm": ("Romanian (small)", "ro", "sm"),
    "ro_core_news_md": ("Romanian (medium)", "ro", "md"),
    "ro_core_news_lg": ("Romanian (large)", "ro", "lg"),

    # Catalan
    "ca_core_news_sm": ("Catalan (small)", "ca", "sm"),
    "ca_core_news_md": ("Catalan (medium)", "ca", "md"),
    "ca_core_news_lg": ("Catalan (large)", "ca", "lg"),

    # Finnish
    "fi_core_news_sm": ("Finnish (small)", "fi", "sm"),
    "fi_core_news_md": ("Finnish (medium)", "fi", "md"),
    "fi_core_news_lg": ("Finnish (large)", "fi", "lg"),

    # Croatian
    "hr_core_news_sm": ("Croatian (small)", "hr", "sm"),
    "hr_core_news_md": ("Croatian (medium)", "hr", "md"),
    "hr_core_news_lg": ("Croatian (large)", "hr", "lg"),

    # Korean
    "ko_core_news_sm": ("Korean (small)", "ko", "sm"),
    "ko_core_news_md": ("Korean (medium)", "ko", "md"),
    "ko_core_news_lg": ("Korean (large)", "ko", "lg"),

    # Macedonian
    "mk_core_news_sm": ("Macedonian (small)", "mk", "sm"),
    "mk_core_news_md": ("Macedonian (medium)", "mk", "md"),
    "mk_core_news_lg": ("Macedonian (large)", "mk", "lg"),

    # Russian
    "ru_core_news_sm": ("Russian (small)", "ru", "sm"),
    "ru_core_news_md": ("Russian (medium)", "ru", "md"),
    "ru_core_news_lg": ("Russian (large)", "ru", "lg"),

    # Ukrainian
    "uk_core_news_sm": ("Ukrainian (small)", "uk", "sm"),
    "uk_core_news_md": ("Ukrainian (medium)", "uk", "md"),
    "uk_core_news_lg": ("Ukrainian (large)", "uk", "lg"),
}


def load_existing_models(data_table_path):
    """Load existing model entries from the data table to avoid duplicates."""
    existing = set()
    if data_table_path and Path(data_table_path).exists():
        with open(data_table_path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if parts:
                        existing.add(parts[0])  # Add the value (first column)
    return existing


def get_model_version(model_name):
    """
    Look up the compatible model version from spaCy's compatibility.json.
    """
    print(f"Looking up compatible version for {model_name} (spaCy {SPACY_VERSION})")
    try:
        response = urllib.request.urlopen(SPACY_COMPAT_URL)
        compat = json.loads(response.read())
        spacy_compat = compat.get("spacy", {})
        if SPACY_VERSION in spacy_compat:
            models = spacy_compat[SPACY_VERSION]
            if model_name in models:
                versions = models[model_name]
                version = versions[0] if isinstance(versions, list) else versions
                print(f"Found compatible version: {version}")
                return version
        print(f"No compatible version found for {model_name} with spaCy {SPACY_VERSION}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error fetching compatibility info: {e}", file=sys.stderr)
        return None


def download_model(model_name):
    """
    Download a spaCy model from GitHub releases using pip install.

    spaCy models are not on PyPI — they are hosted on GitHub releases.
    This function looks up the compatible version and installs from the
    direct URL, without requiring spaCy itself.
    """
    print(f"Downloading spaCy model: {model_name}")

    version = get_model_version(model_name)
    if not version:
        return False

    # Construct GitHub releases URL
    wheel_name = f"{model_name}-{version}-py3-none-any.whl"
    url = f"{SPACY_MODELS_BASE_URL}/{model_name}-{version}/{wheel_name}"
    print(f"Installing from {url}")

    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", url],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        print(f"Successfully downloaded {model_name} {version}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error downloading model: {e}", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Download and register spaCy language models")
    parser.add_argument("--model", action="append", required=True, choices=SPACY_MODELS.keys(),
                        help="spaCy model(s) to download (can be specified multiple times)")
    parser.add_argument("--output", required=True,
                        help="JSON output file for Galaxy data manager")
    parser.add_argument("--data-table", required=False,
                        help="Path to existing data table file to check for duplicates")

    args = parser.parse_args()

    # Load existing models to avoid duplicates
    existing_models = load_existing_models(args.data_table)

    # List to collect all data table entries
    data_table_entries = []

    # Process each model
    for model_name in args.model:
        if model_name in existing_models:
            print(f"\n{'=' * 60}")
            print(f"Skipping {model_name} - already in data table")
            print(f"{'=' * 60}")
        else:
            print(f"\n{'=' * 60}")
            print(f"Processing {model_name}...")
            print(f"{'=' * 60}")

            display_name, language, size = SPACY_MODELS[model_name]

            # Download model via pip to verify it exists and is installable
            if not download_model(model_name):
                print(f"WARNING: Failed to download {model_name}", file=sys.stderr)
                continue  # Skip this model but continue with others

            # Verify the package was installed
            try:
                subprocess.run(
                    [sys.executable, "-m", "pip", "show", model_name],
                    capture_output=True,
                    text=True,
                    check=True
                )
                print(f"Model package {model_name} verified")
            except subprocess.CalledProcessError:
                print(f"WARNING: Package {model_name} not found after install", file=sys.stderr)
                continue

            # Add to data table entries
            data_table_entries.append({
                "value": model_name,
                "name": display_name,
                "lang": language,
                "size": size,
                "model_name": model_name
            })

            print(f"Successfully registered {display_name}")
            print(f"  Model name: {model_name}")
            print(f"  Language: {language}")
            print(f"  Size: {size}")

    # Create data manager JSON output
    data_manager_output = {
        "data_tables": {
            "spacy_models": data_table_entries
        }
    }

    # Write output JSON
    with open(args.output, "w") as f:
        json.dump(data_manager_output, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"Summary: Successfully registered {len(data_table_entries)} model(s)")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
