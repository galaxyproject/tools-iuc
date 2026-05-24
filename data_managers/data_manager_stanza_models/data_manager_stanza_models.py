#!/usr/bin/env python
"""
Data Manager for Stanza Language Models

Downloads Stanza language models from HuggingFace and registers them in
Galaxy's data table. Does NOT require stanza to be installed — downloads
the default model package (zip) directly via HTTP.
"""

import argparse
import json
import sys
import urllib.request
import zipfile
from pathlib import Path


# Stanza resource configuration
STANZA_VERSION = "1.11.0"
RESOURCES_URL = f"https://raw.githubusercontent.com/stanfordnlp/stanza-resources/main/resources_{STANZA_VERSION}.json"
# URL template: filled with lang and resources_version from the resources JSON
DEFAULT_URL_TEMPLATE = "https://huggingface.co/stanfordnlp/stanza-{lang}/resolve/v{resources_version}/models/{filename}"


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


def fetch_resources():
    """Fetch the Stanza resources JSON to get download URLs and checksums."""
    print(f"Fetching Stanza resources from {RESOURCES_URL}")
    response = urllib.request.urlopen(RESOURCES_URL)
    return json.loads(response.read())


def download_model(lang, model_dir, resources):
    """Download a Stanza language model package from HuggingFace.

    Downloads the default.zip package for the language and extracts it
    into the model_dir/<lang>/ directory. Also writes the resources.json
    file needed by Stanza at runtime.
    """
    # Get the URL template from the resources JSON
    url_template = resources.get("url", DEFAULT_URL_TEMPLATE)

    # Check if the language exists in resources
    if lang not in resources:
        print(f"Error: Language '{lang}' not found in Stanza resources", file=sys.stderr)
        return False

    # Download the default_fast.zip package (nocharlm models — much lower memory usage)
    # Fall back to default.zip if default_fast is not available for this language
    packages = resources.get(lang, {}).get("packages", {})
    package_name = "default_fast" if "default_fast" in packages else "default"
    zip_url = url_template.format(
        lang=lang,
        resources_version=STANZA_VERSION,
        filename=f"{package_name}.zip"
    )
    print(f"Using package: {package_name}")

    lang_dir = Path(model_dir) / lang
    lang_dir.mkdir(parents=True, exist_ok=True)
    zip_path = lang_dir / "default.zip"

    print(f"Downloading {zip_url}")
    try:
        urllib.request.urlretrieve(zip_url, str(zip_path))
    except Exception as e:
        print(f"Error downloading {lang} model: {e}", file=sys.stderr)
        return False

    # Extract the zip
    print(f"Extracting to {lang_dir}")
    try:
        with zipfile.ZipFile(str(zip_path), 'r') as zf:
            zf.extractall(str(lang_dir))
    except Exception as e:
        print(f"Error extracting {lang} model: {e}", file=sys.stderr)
        return False

    # Write resources.json if it doesn't exist yet (needed by stanza.Pipeline)
    resources_path = Path(model_dir) / "resources.json"
    if resources_path.exists():
        with open(resources_path) as f:
            existing_resources = json.load(f)
    else:
        existing_resources = {}

    # Add/update this language's resource entry
    existing_resources[lang] = resources[lang]
    # Also include the URL key
    existing_resources["url"] = url_template

    with open(resources_path, 'w') as f:
        json.dump(existing_resources, f, indent=2)

    # Clean up the zip file
    zip_path.unlink()

    print(f"Successfully downloaded and extracted {lang} model")
    return True


def main():
    parser = argparse.ArgumentParser(description="Download and register Stanza language models")
    parser.add_argument("--model", action="append", required=True,
                        help="Language code(s) to download (can be specified multiple times)")
    parser.add_argument("--target-directory", required=True,
                        help="Persistent directory to store downloaded models")
    parser.add_argument("--output", required=True,
                        help="JSON output file for Galaxy data manager")
    parser.add_argument("--data-table", required=False,
                        help="Path to existing data table file to check for duplicates")

    args = parser.parse_args()

    # Load existing models to avoid duplicates
    existing_models = load_existing_models(args.data_table)

    # Fetch resources JSON
    try:
        resources = fetch_resources()
    except Exception as e:
        print(f"Error fetching Stanza resources: {e}", file=sys.stderr)
        sys.exit(1)

    # Use the persistent target directory for models
    model_dir = Path(args.target_directory)
    model_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    for lang in args.model:
        if lang in existing_models:
            print(f"\n{'=' * 60}")
            print(f"Skipping {lang} - already in data table")
            print(f"{'=' * 60}")
            continue

        print(f"\n{'=' * 60}")
        print(f"Processing {lang}...")
        print(f"{'=' * 60}")

        display_name = STANZA_LANGUAGES.get(lang, lang)

        if not download_model(lang, model_dir, resources):
            print(f"WARNING: Failed to download {lang}", file=sys.stderr)
            continue

        data_table_entries.append({
            "value": lang,
            "name": display_name,
            "lang": lang,
            "models_path": str(model_dir),
        })

        print(f"Successfully registered {display_name}")
        print(f"  Language code: {lang}")
        print(f"  Models path: {model_dir}")

    # Create data manager JSON output
    data_manager_output = {
        "data_tables": {
            "stanza_models": data_table_entries
        }
    }

    with open(args.output, "w") as f:
        json.dump(data_manager_output, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"Summary: Successfully registered {len(data_table_entries)} model(s)")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
