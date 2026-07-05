#!/usr/bin/env python
"""
Data Manager for Stanza Language Models

Downloads Stanza language models with ``stanza.download()`` into the data
manager output's extra files directory and registers them in Galaxy's
``stanza_models`` data table. Galaxy then moves the models into the managed
data directory (see data_manager_conf.xml).

Each language is installed as a self-contained stanza_resources directory
(resources.json plus a per-language model subdirectory) using the
memory-efficient ``default_fast`` package, matching the Stanza NLP tool which
loads it via ``stanza.Pipeline(dir=models_path, package="default_fast")``.
"""

import argparse
import json
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


def main():
    parser = argparse.ArgumentParser(description="Download and register Stanza language models")
    parser.add_argument("data_manager_json",
                        help="Galaxy data manager JSON file (prepopulated with output paths)")
    parser.add_argument("--model", action="append", required=True,
                        help="Language code(s) to download (can be specified multiple times)")
    args = parser.parse_args()

    # Galaxy prepopulates the data manager JSON with the output dataset's
    # extra files path, which is a writable directory. Models are downloaded
    # there and Galaxy moves them into the managed data directory afterwards.
    with open(args.data_manager_json) as fh:
        params = json.load(fh)
    target_dir = Path(params["output_data"][0]["extra_files_path"])
    target_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    for lang in args.model:
        display_name = STANZA_LANGUAGES.get(lang, lang)

        # Install each language as its own stanza_resources directory so the
        # tool can load it directly with stanza.Pipeline(dir=models_path).
        lang_dir = target_dir / lang
        print(f"Downloading {display_name} ({lang}) models with the default_fast package...")
        stanza.download(
            lang=lang,
            model_dir=str(lang_dir),
            package="default_fast",
            verbose=False,
        )

        data_table_entries.append({
            "value": lang,
            "name": display_name,
            "lang": lang,
            # Relative to extra_files_path; data_manager_conf.xml moves it into
            # ${GALAXY_DATA_MANAGER_DATA_PATH}/stanza_models/${value}.
            "models_path": lang,
        })

        print(f"Successfully downloaded {display_name} ({lang})")

    data_manager_output = {
        "data_tables": {
            "stanza_models": data_table_entries
        }
    }

    with open(args.data_manager_json, "w") as fh:
        json.dump(data_manager_output, fh, sort_keys=True)

    print(f"Summary: Successfully registered {len(data_table_entries)} model(s)")


if __name__ == "__main__":
    main()
