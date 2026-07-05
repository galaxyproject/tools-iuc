#!/usr/bin/env python
"""
Data Manager for Stanza Language Models

Downloads Stanza language models with ``stanza.download()`` and registers them
in Galaxy's ``stanza_models`` data table. Uses the memory-efficient
``default_fast`` package (nocharlm models), matching the Stanza NLP tool so the
same on-disk layout (a stanza_resources directory with resources.json) is loaded
at run time via ``stanza.Pipeline(dir=models_path, package="default_fast")``.
"""

import argparse
import json
import sys
from pathlib import Path

import stanza


# Language display names
STANZA_LANGUAGES = {
    "en": "English",
    "zh-hans": "Chinese (Simplified)",
    "zh-hant": "Chinese (Traditional)",
    "ar": "Arabic",
    "fr": "French",
    "de": "German",
    "es": "Spanish",
    "it": "Italian",
    "pt": "Portuguese",
    "nl": "Dutch",
    "ru": "Russian",
    "uk": "Ukrainian",
    "pl": "Polish",
    "ja": "Japanese",
    "ko": "Korean",
    "hi": "Hindi",
    "tr": "Turkish",
    "el": "Greek",
    "hu": "Hungarian",
    "sv": "Swedish",
    "da": "Danish",
    "nb": "Norwegian Bokmål",
    "nn": "Norwegian Nynorsk",
    "fi": "Finnish",
    "ro": "Romanian",
    "ca": "Catalan",
    "cs": "Czech",
    "sk": "Slovak",
    "sl": "Slovenian",
    "hr": "Croatian",
    "sr": "Serbian",
    "bg": "Bulgarian",
    "lv": "Latvian",
    "lt": "Lithuanian",
    "et": "Estonian",
    "he": "Hebrew",
    "fa": "Persian",
    "vi": "Vietnamese",
    "th": "Thai",
    "id": "Indonesian",
    "af": "Afrikaans",
    "eu": "Basque",
    "gl": "Galician",
    "hy": "Armenian",
    "ka": "Georgian",
    "ta": "Tamil",
    "te": "Telugu",
    "mr": "Marathi",
    "ur": "Urdu",
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
                        existing.add(parts[0])
    return existing


def main():
    parser = argparse.ArgumentParser(description="Download and register Stanza language models")
    parser.add_argument("--model", action="append", required=True,
                        help="Language code(s) to download (can be specified multiple times)")
    parser.add_argument("--target-directory", required=True,
                        help="Persistent stanza_resources directory to store downloaded models")
    parser.add_argument("--output", required=True,
                        help="JSON output file for Galaxy data manager")
    parser.add_argument("--data-table", required=False,
                        help="Path to existing data table file to check for duplicates")

    args = parser.parse_args()

    # Load existing models to avoid duplicates
    existing_models = load_existing_models(args.data_table)

    # Use the persistent target directory as the stanza_resources root. stanza
    # writes resources.json here plus a per-language subdirectory.
    model_dir = Path(args.target_directory)
    model_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    for lang in args.model:
        if lang in existing_models:
            print(f"Skipping {lang} - already in data table")
            continue

        display_name = STANZA_LANGUAGES.get(lang, lang)

        print(f"Downloading {display_name} ({lang}) models with the default_fast package...")
        try:
            stanza.download(
                lang=lang,
                model_dir=str(model_dir),
                package="default_fast",
                verbose=False,
            )
        except Exception as e:
            print(f"Error downloading {lang} model: {e}", file=sys.stderr)
            sys.exit(1)

        data_table_entries.append({
            "value": lang,
            "name": display_name,
            "lang": lang,
            "models_path": str(model_dir),
        })

        print(f"Successfully registered {display_name} ({lang}) at {model_dir}")

    # Create data manager JSON output
    data_manager_output = {
        "data_tables": {
            "stanza_models": data_table_entries
        }
    }

    with open(args.output, "w") as f:
        json.dump(data_manager_output, f, indent=2)

    print(f"Summary: Successfully registered {len(data_table_entries)} model(s)")


if __name__ == "__main__":
    main()
