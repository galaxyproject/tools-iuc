#!/usr/bin/env python
"""
Data Manager for Stanza Language Models

Downloads Stanza language models into the data manager output's extra files
directory and registers them in Galaxy's ``stanza_models`` data table. Galaxy
then moves the models into the managed data directory (see data_manager_conf.xml).

The ``default`` package is downloaded with ``stanza.download()`` (single zip
archive). ``default_fast`` and ``default_accurate`` are downloaded file-by-file
with stanza's own package-expansion helpers, because stanza.download() crashes
on languages whose packages use the identity lemmatizer (see
download_language_models).

Each language is installed as a self-contained stanza_resources directory
(resources.json plus a per-language model subdirectory) using the selected
package (``default_fast`` by default). The Stanza NLP tool loads it via
``stanza.Pipeline(dir=models_path, package=<package>)``, reading the package
back from the data table. Installing the ``default`` or ``default_accurate``
package additionally provides constituency parsing, which ``default_fast``
does not include.
"""

import argparse
import json
import os
import shutil
from pathlib import Path

import stanza
from certifi import where as certifi_where
from stanza._version import __resources_version__, __version__
from stanza.resources.common import (
    add_dependencies,
    download_resources_json,
    flatten_processor_list,
    LEMMA,
    load_resources_json,
    maintain_processor_list,
    PROCESSOR_VARIANTS,
    request_file,
)

# Ensure CA certificates are found by hf_xet in conda containers.
os.environ.setdefault("SSL_CERT_FILE", certifi_where())


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


def _resolve_alias(resources, lang):
    """Follow 'alias' entries in the resources file (mirrors stanza.download)."""
    while isinstance(resources.get(lang), dict) and 'alias' in resources[lang]:
        lang = resources[lang]['alias']
    return lang


def download_language_models(lang, model_dir, package, resources_version, resources):
    """
    Download the models of one language for a package that is distributed as
    individual model files (default_fast, default_accurate).

    stanza.download() cannot be used for these packages: languages whose
    packages use the identity lemmatizer (lemma == token; currently Vietnamese,
    Thai, Telugu, Maltese and Faroese) make its download_models() look up
    resources[lang]['lemma']['identity']['md5'], but the identity lemmatizer
    ships no model file and has no entry in resources.json. Stanza therefore
    raises "Cannot find the following processor and model name combination:
    lemma, identity" (stanza >= 1.12 through current main). The 'default'
    package is unaffected because it downloads a single default.zip instead.

    This replicates the per-file download using stanza's own package-expansion
    helpers and checksum-verified request_file, skipping processors that are
    implemented in code and need no model file. The pipeline side handles these
    variants itself (e.g. the lemma_use_identity option), so the resulting
    stanza_resources directory loads exactly like a stanza.download() one.
    """
    lang = _resolve_alias(resources, lang)
    download_list = maintain_processor_list(resources, lang, package, {}, allow_pretrain=True)
    download_list = add_dependencies(resources, lang, download_list)
    download_list = flatten_processor_list(download_list)

    url = resources['url']
    for processor, model_name in download_list:
        # External processor variants (spacy, jieba, sudachipy, ...) and the
        # identity lemmatizer are code-only and ship no downloadable model.
        if model_name in PROCESSOR_VARIANTS.get(processor, {}):
            continue
        if processor == LEMMA and model_name == 'identity':
            continue
        request_file(
            url.format(resources_version=resources_version, lang=lang,
                       filename=f"{processor}/{model_name}.pt"),
            os.path.join(model_dir, lang, processor, f'{model_name}.pt'),
            None,
            md5=resources[lang][processor][model_name]['md5'],
            alternate_md5=resources[lang][processor][model_name].get('alternate_md5', None),
        )


def main():
    parser = argparse.ArgumentParser(description="Download and register Stanza language models")
    parser.add_argument("data_manager_json",
                        help="Galaxy data manager JSON file (prepopulated with output paths)")
    parser.add_argument("--model", action="append", required=True,
                        help="Language code(s) to download (can be specified multiple times)")
    parser.add_argument("--package", default="default_fast",
                        help="Stanza model package to install (default_fast, default, or "
                             "default_accurate). default/default_accurate add constituency "
                             "parsing at the cost of larger downloads.")
    parser.add_argument("--stanza-version", default=__version__,
                        help="Stanza library version the models are downloaded with "
                             "(e.g. 1.12.0). Defaults to the installed stanza library "
                             "version.")
    parser.add_argument("--resources-version", default=__resources_version__,
                        help="Stanza model resources version to download (e.g. 1.12.0). "
                             "Defaults to the resources version of the installed stanza "
                             "library, which is distinct from the library version and may "
                             "diverge from it.")
    args = parser.parse_args()

    # Galaxy prepopulates the data manager JSON with the output dataset's
    # extra files path, which is a writable directory. Models are downloaded
    # there and Galaxy moves them into the managed data directory afterwards.
    with open(args.data_manager_json) as fh:
        params = json.load(fh)
    target_dir = Path(params["output_data"][0]["extra_files_path"])
    target_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    package = args.package
    stanza_version = args.stanza_version
    resources_version = args.resources_version

    # The package expansions and model checksums below come from resources.json,
    # so fetch it once up front and reuse it for every selected language. It is
    # copied into each per-language directory, which the Stanza NLP tool loads
    # as a self-contained stanza_resources directory.
    download_resources_json(str(target_dir), resources_version=resources_version)
    resources = load_resources_json(str(target_dir))
    resources_json = target_dir / "resources.json"

    for lang in args.model:
        display_name = STANZA_LANGUAGES.get(lang, lang)

        # A data-table entry is keyed by language + package + resources version so
        # multiple packages and versions for the same language can coexist.
        # The value doubles as the on-disk subdirectory name.
        entry_value = f"{lang}-{package}-v{resources_version}"
        lang_dir = target_dir / entry_value
        lang_dir.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {display_name} ({lang}) models with the {package} package...")

        if package == "default":
            # The 'default' package is distributed as a single zip archive;
            # stanza.download() handles it (and the resources.json copy) itself.
            stanza.download(
                lang=lang,
                model_dir=str(lang_dir),
                package=package,
                resources_version=resources_version,
                verbose=False,
            )
        else:
            # Per-file packages: stanza.download() crashes on identity-lemma
            # languages (see download_language_models), so download the model
            # files ourselves.
            shutil.copy(resources_json, lang_dir / "resources.json")
            download_language_models(lang, str(lang_dir), package, resources_version, resources)

        data_table_entries.append({
            "value": entry_value,
            "name": f"{display_name} — {package}",
            "lang": lang,
            "package": package,
            "stanza_version": stanza_version,
            "resources_version": resources_version,
            # Relative to extra_files_path; data_manager_conf.xml moves it into
            # ${GALAXY_DATA_MANAGER_DATA_PATH}/stanza_models/${value}.
            "models_path": entry_value,
        })

        print(f"Successfully downloaded {display_name} ({lang}) [{package}]")

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
