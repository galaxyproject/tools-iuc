#!/usr/bin/env python
"""
Data Manager for Stanford CoreNLP Language Models

Downloads CoreNLP model JARs from Maven Central into the data manager output's
extra files path and registers them in Galaxy's ``corenlp_models`` data table.
Galaxy then moves the JARs into a shared managed directory (see
data_manager_conf.xml) so the CoreNLP tool can symlink them at run time. All
JARs land in the same directory so the coreference annotator can find the
common models JAR alongside a language JAR.
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path

CORENLP_VERSION = "4.5.10"
BASE_URL = f"https://repo1.maven.org/maven2/edu/stanford/nlp/stanford-corenlp/{CORENLP_VERSION}"

# Common models JAR (dictionaries/models needed for coreference).
COMMON_MODELS = {
    "name": "Common Models",
    "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models.jar",
    "url": f"{BASE_URL}/stanford-corenlp-{CORENLP_VERSION}-models.jar",
}

LANGUAGE_MODELS = {
    "ar": {"name": "Arabic", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-arabic.jar"},
    "zh": {"name": "Chinese", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-chinese.jar"},
    "en": {"name": "English", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-english.jar"},
    "fr": {"name": "French", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-french.jar"},
    "de": {"name": "German", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-german.jar"},
    "hu": {"name": "Hungarian", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-hungarian.jar"},
    "it": {"name": "Italian", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-italian.jar"},
    "es": {"name": "Spanish", "jar_name": f"stanford-corenlp-{CORENLP_VERSION}-models-spanish.jar"},
}
for _code, _info in LANGUAGE_MODELS.items():
    _info["url"] = f"{BASE_URL}/{_info['jar_name']}"

# Lightweight stand-in used only by the data manager test (--test): the ~8 MB
# CoreNLP code JAR, so the test exercises download and registration without
# pulling a multi-hundred-MB model JAR.
TEST_MODEL = {
    "value": "en",
    "name": "English",
    "lang_code": "en",
    "jar_name": f"stanford-corenlp-{CORENLP_VERSION}.jar",
    "url": f"{BASE_URL}/stanford-corenlp-{CORENLP_VERSION}.jar",
}


def download_model(url, target_path):
    """Download a file from URL to target_path."""
    print(f"Downloading {url}")
    print(f"Saving to {target_path}")
    try:
        urllib.request.urlretrieve(url, target_path)
        print("Download complete")
        return True
    except Exception as e:
        print(f"Error downloading file: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Download and register CoreNLP language models")
    parser.add_argument("data_manager_json",
                        help="Galaxy data manager JSON file (prepopulated with output paths)")
    parser.add_argument("--language", action="append", choices=LANGUAGE_MODELS.keys(),
                        help="Language code(s) for the model(s) to download")
    parser.add_argument("--common-models", action="store_true",
                        help="Download the common models JAR (required for coreference)")
    parser.add_argument("--test", action="store_true",
                        help="Download the small CoreNLP code JAR as a lightweight test stand-in")
    args = parser.parse_args()

    if not args.language and not args.common_models and not args.test:
        parser.error("one of --language, --common-models or --test is required")

    # Galaxy prepopulates the data manager JSON with the output dataset's extra
    # files path, a writable directory. JARs are downloaded there and Galaxy
    # moves them into the managed data directory afterwards.
    with open(args.data_manager_json) as fh:
        params = json.load(fh)
    target_dir = Path(params["output_data"][0]["extra_files_path"])
    target_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    def install(value, name, lang_code, jar_name, url):
        if not download_model(url, str(target_dir / jar_name)):
            return False
        data_table_entries.append({
            "value": value,
            "name": name,
            "lang_code": lang_code,
            "version": CORENLP_VERSION,
            # Relative to extra_files_path; data_manager_conf.xml moves it into
            # the shared ${GALAXY_DATA_MANAGER_DATA_PATH}/corenlp/<version> directory.
            "models_path": jar_name,
        })
        print(f"Registered {name}")
        return True

    if args.test:
        if not install(TEST_MODEL["value"], TEST_MODEL["name"], TEST_MODEL["lang_code"],
                       TEST_MODEL["jar_name"], TEST_MODEL["url"]):
            sys.exit(1)
    else:
        if args.common_models:
            if not install("common", COMMON_MODELS["name"], "common",
                           COMMON_MODELS["jar_name"], COMMON_MODELS["url"]):
                sys.exit(1)
        for lang_code in args.language or []:
            info = LANGUAGE_MODELS[lang_code]
            if not install(lang_code, info["name"], lang_code, info["jar_name"], info["url"]):
                sys.exit(1)

    data_manager_output = {"data_tables": {"corenlp_models": data_table_entries}}
    with open(args.data_manager_json, "w") as fh:
        json.dump(data_manager_output, fh, sort_keys=True)

    print(f"Summary: {len(data_table_entries)} model(s) registered")


if __name__ == "__main__":
    main()
