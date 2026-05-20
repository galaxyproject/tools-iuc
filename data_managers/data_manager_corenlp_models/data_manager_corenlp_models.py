#!/usr/bin/env python
# Copyright 2006 The Galaxy Project. All rights reserved.
"""
Data Manager for Stanford CoreNLP Language Models

Downloads CoreNLP language model JARs to a persistent directory and registers
them in the Galaxy data table. JARs are stored at the absolute path so the
CoreNLP tool can symlink them at runtime.
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path


# CoreNLP version and model information
CORENLP_VERSION = "4.5.10"

# Common models JAR (contains dcoref dictionaries and common models)
COMMON_MODELS = {
    "name": "Common Models",
    "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models.jar",
    "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models.jar"
}

LANGUAGE_MODELS = {
    "ar": {
        "name": "Arabic",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-arabic.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-arabic.jar"
    },
    "zh": {
        "name": "Chinese",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-chinese.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-chinese.jar"
    },
    "en": {
        "name": "English",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-english.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-english.jar"
    },
    "fr": {
        "name": "French",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-french.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-french.jar"
    },
    "de": {
        "name": "German",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-german.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-german.jar"
    },
    "hu": {
        "name": "Hungarian",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-hungarian.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-hungarian.jar"
    },
    "it": {
        "name": "Italian",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-italian.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-italian.jar"
    },
    "es": {
        "name": "Spanish",
        "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-spanish.jar",
        "url": f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}/stanford-corenlp-{CORENLP_VERSION}-models-spanish.jar"
    }
}


def download_model(url, target_path):
    """Download a file from URL to target path with progress reporting."""
    print(f"Downloading from {url}")
    print(f"Saving to {target_path}")

    def report_progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(100, (downloaded / total_size) * 100)
            mb = downloaded / 1024 / 1024
            total_mb = total_size / 1024 / 1024
            print(f"  {percent:.0f}% ({mb:.0f}/{total_mb:.0f} MB)", flush=True)

    try:
        urllib.request.urlretrieve(url, target_path, reporthook=report_progress)
        print("Download complete!")
        return True
    except Exception as e:
        print(f"Error downloading file: {e}", file=sys.stderr)
        return False


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
                        existing.add(parts[0])
    return existing


def main():
    parser = argparse.ArgumentParser(description="Download and register CoreNLP language models")
    parser.add_argument("--language", action="append", choices=LANGUAGE_MODELS.keys(),
                        help="Language code(s) for the model(s) to download")
    parser.add_argument("--common-models", action="store_true",
                        help="Download the common models JAR (required for coreference)")
    parser.add_argument("--target-directory", required=True,
                        help="Directory to store the downloaded model JARs")
    parser.add_argument("--output", required=True,
                        help="JSON output file for Galaxy data manager")
    parser.add_argument("--data-table", required=False,
                        help="Path to existing data table file to check for duplicates")

    args = parser.parse_args()

    if not args.language and not args.common_models:
        parser.error("At least one of --language or --common-models must be specified")

    existing_models = load_existing_models(args.data_table)

    target_dir = Path(args.target_directory)
    target_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    # Process common models if requested
    if args.common_models:
        if "common" in existing_models:
            print(f"\n{'=' * 60}")
            print(f"Skipping {COMMON_MODELS['name']} - already in data table")
            print(f"{'=' * 60}")
        else:
            print(f"\n{'=' * 60}")
            print(f"Processing {COMMON_MODELS['name']}...")
            print(f"{'=' * 60}")

            jar_path = target_dir / COMMON_MODELS["jar_name"]

            if jar_path.exists():
                print(f"Model already exists at {jar_path}")
            else:
                if not download_model(COMMON_MODELS["url"], str(jar_path)):
                    print(f"WARNING: Failed to download {COMMON_MODELS['name']}", file=sys.stderr)

            if jar_path.exists():
                data_table_entries.append({
                    "value": "common",
                    "name": COMMON_MODELS["name"],
                    "lang_code": "common",
                    "models_path": str(jar_path.absolute())
                })
                print(f"Registered {COMMON_MODELS['name']}")
                print(f"  Path: {jar_path.absolute()}")

    # Process each language
    if args.language:
        for lang_code in args.language:
            if lang_code in existing_models:
                print(f"\n{'=' * 60}")
                print(f"Skipping {LANGUAGE_MODELS[lang_code]['name']} - already in data table")
                print(f"{'=' * 60}")
            else:
                model_info = LANGUAGE_MODELS[lang_code]

                print(f"\n{'=' * 60}")
                print(f"Processing {model_info['name']} model...")
                print(f"{'=' * 60}")

                jar_path = target_dir / model_info["jar_name"]

                if jar_path.exists():
                    print(f"Model already exists at {jar_path}")
                else:
                    if not download_model(model_info["url"], str(jar_path)):
                        print(f"WARNING: Failed to download {model_info['name']} model", file=sys.stderr)
                        continue

                data_table_entries.append({
                    "value": lang_code,
                    "name": model_info["name"],
                    "lang_code": lang_code,
                    "models_path": str(jar_path.absolute())
                })
                print(f"Registered {model_info['name']} model")
                print(f"  Path: {jar_path.absolute()}")

    # Write data manager JSON output
    data_manager_output = {
        "data_tables": {
            "corenlp_models": data_table_entries
        }
    }

    with open(args.output, "w") as f:
        json.dump(data_manager_output, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"Summary: {len(data_table_entries)} model(s) registered")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
