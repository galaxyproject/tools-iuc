#!/usr/bin/env python
"""
Data Manager for spaCy Language Models

Downloads spaCy model wheels from the explosion/spacy-models GitHub releases,
extracts the model data directory into the data manager output's extra files
path, and registers the path in Galaxy's ``spacy_models`` data table. Galaxy
then moves the model into the managed data directory (see data_manager_conf.xml).

Storing the extracted model directory (rather than pip-installing the package)
means the tool can load it directly and offline with ``spacy.load(model_path)``,
so the model persists across jobs and does not need to be re-downloaded at run
time. No spaCy install is required to run this data manager.
"""

import argparse
import json
import shutil
import urllib.request
import zipfile
from pathlib import Path
import sys

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


def get_model_version(model_name):
    """Look up the compatible model version from spaCy's compatibility.json."""
    print(f"Looking up compatible version for {model_name} (spaCy {SPACY_VERSION})")
    response = urllib.request.urlopen(SPACY_COMPAT_URL)
    compat = json.loads(response.read())
    models = compat.get("spacy", {}).get(SPACY_VERSION, {})
    versions = models.get(model_name)
    if not versions:
        raise RuntimeError(f"No compatible version for {model_name} with spaCy {SPACY_VERSION}")
    version = versions[0] if isinstance(versions, list) else versions
    print(f"Found compatible version: {version}")
    return version

def install_wheel(wheel_url, dest_dir):
    """Download a spaCy model wheel and extract its data directory into dest_dir.

    spaCy model wheels contain ``<lang>_<name>/<lang>_<name>-<version>/`` which is
    the loadable model data directory (config.cfg, meta.json, pipeline
    components). Its contents are written to ``dest_dir/<lang>_<name>_<version>``
     so that ``spacy.load(dest_dir/<lang>_<name>_<version>)`` works offline.
    """

    dest_dir = Path(dest_dir)
    wheel_path = dest_dir / "model.whl"
    extract_dir = dest_dir / "extracted_wheel"

    print(f"Downloading {wheel_url}")
    urllib.request.urlretrieve(wheel_url, str(wheel_path))

    print(f"Extracting {wheel_path}")
    with zipfile.ZipFile(str(wheel_path)) as zf:
        zf.extractall(str(extract_dir))

    # find all the meta.json files in the wheel and check
    # that there is only 1 --
    meta_file = list(extract_dir.glob("*/meta.json"))
    if len(meta_file) == 0:
        raise RuntimeError("Couldn't find meta.json. Is this really a spaCy model?")
    elif len(meta_file) != 1:
        raise RuntimeError("Found multiple meta.json mfiles. Is this really a spaCy model?")

    # read the metadata out of the file
    with open(meta_file[0]) as f:
        metadata = json.load(f)

    # quick sanity check to ensure that this metadata really
    # is describing a spaCy model.
    # NOTE: ideally we should also check the version range
    #       against the SPACY_VERSION constant, but as the
    #       user is manually installing this version we are
    #       currently trusting them to know what they are doing
    if "spacy_version" not in metadata:
        raise RuntimeError("spacy_version not specified in meta.json. Is this wheel really a spaCy model?")

    # get the lang, name, and version from the metadat and
    # build the appropriate variables we want later
    model_name = f"{metadata["lang"]}_{metadata["name"]}"
    version = metadata["version"]

    print(f"Extracting {model_name} v{version}")

    # this is the data directory inside the extracted wheel, if it
    # doesn't exist then we have a big problem
    data_dir = extract_dir / model_name / f"{model_name}-{version}"
    if not data_dir.is_dir():
        raise RuntimeError(f"Model data directory not found in wheel: {data_dir}")

    # move the extracted data directory to where Galaxy can find it
    dest_dir = dest_dir / f"{model_name}_{version}"
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    shutil.move(str(data_dir), str(dest_dir))

    # Clean up the wheel and leftover extraction directory.
    wheel_path.unlink()
    shutil.rmtree(extract_dir, ignore_errors=True)

    # I don't like this but... for consistency if we have installed a wheel for
    # one of the standard models then use the same name so that we avoid having the
    # same model installed twice under different descriptions
    desc = SPACY_MODELS[model_name][0] if model_name in SPACY_MODELS else model_name

    print(f"Successfully installed {desc} v{version}")
    
    return {
        "value": f"{model_name}_{version}",
        "name": f"{desc} v{version}",
        "model": metadata["name"],
        "lang": metadata["lang"],
        "version": version,
        # Relative to extra_files_path; data_manager_conf.xml moves it into
        # ${GALAXY_DATA_MANAGER_DATA_PATH}/spacy_models/${value}/
        "model_path": f"{model_name}_{version}"
    }
    

def install_model(model_name, dest_dir):
    """Download a standard spaCy model wheel and extract its data directory into dest_dir.

    The exact version of the model to download is determined by checking the
    compatiability list using the language and model name

    spaCy model wheels contain ``<lang>_<name>/<lang>_<name>-<version>/`` which is
    the loadable model data directory (config.cfg, meta.json, pipeline
    components). Its contents are written to ``dest_dir/<lang>_<name>_<version>``
     so that ``spacy.load(dest_dir/<lang>_<name>_<version>)`` works offline.
    """
    version = get_model_version(model_name)
    wheel_name = f"{model_name}-{version}-py3-none-any.whl"
    url = f"{SPACY_MODELS_BASE_URL}/{model_name}-{version}/{wheel_name}"

    dest_dir = Path(dest_dir)
    extract_dir = dest_dir / f"_{model_name}_wheel"
    wheel_path = dest_dir / wheel_name

    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, str(wheel_path))

    print(f"Extracting {wheel_name}")
    with zipfile.ZipFile(str(wheel_path)) as zf:
        zf.extractall(str(extract_dir))

    # The loadable model data directory is <model_name>/<model_name>-<version>.
    data_dir = extract_dir / model_name / f"{model_name}-{version}"
    if not data_dir.is_dir():
        raise RuntimeError(f"Model data directory not found in wheel: {data_dir}")

    dest_dir = dest_dir / f"{model_name}_{version}"

    # Move the model data directory contents to dest_dir so dest_dir is the
    # directory that spacy.load() consumes.
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    shutil.move(str(data_dir), str(dest_dir))

    # Clean up the wheel and leftover extraction directory.
    wheel_path.unlink()
    shutil.rmtree(extract_dir, ignore_errors=True)

    return version


def main():
    parser = argparse.ArgumentParser(description="Download and register spaCy language models")
    parser.add_argument("data_manager_json",
                        help="Galaxy data manager JSON file (prepopulated with output paths)")
    
    parser.add_argument("--model", action="append", choices=SPACY_MODELS.keys(),
                        help="spaCy model(s) to download (can be specified multiple times)")
    parser.add_argument("--wheel", action="append", help="spaCy model wheel(s) to download (can be specified multiple times)")
    
    args = parser.parse_args()

    # Galaxy prepopulates the data manager JSON with the output dataset's extra
    # files path, a writable directory. Models are extracted there and Galaxy
    # moves them into the managed data directory afterwards.
    with open(args.data_manager_json) as fh:
        params = json.load(fh)
    target_dir = Path(params["output_data"][0]["extra_files_path"])
    target_dir.mkdir(parents=True, exist_ok=True)

    data_table_entries = []

    if args.model:
        for model_name in args.model:
            display_name, language, size = SPACY_MODELS[model_name]
            print(f"\n{'=' * 60}\nProcessing {model_name}...\n{'=' * 60}")

            version = install_model(model_name, target_dir)

            data_table_entries.append({
                "value": f"{model_name}_{version}",
                "name": f"{display_name} v{version}",
                "model": model_name[3:],
                "lang": language,
                "version": version,
                # Relative to extra_files_path; data_manager_conf.xml moves it into
                # ${GALAXY_DATA_MANAGER_DATA_PATH}/spacy_models/${value}.
                "model_path": f"{model_name}_{version}",
            })

            print(f"Successfully installed {display_name} v{version}")
    
    if args.wheel:
        for wheel_url in args.wheel:
            print(f"\n{'=' * 60}\nProcessing...\n{'=' * 60}")
            data_table_entries.append(install_wheel(wheel_url, target_dir))

    data_manager_output = {"data_tables": {"spacy_models": data_table_entries}}
    with open(args.data_manager_json, "w") as fh:
        json.dump(data_manager_output, fh, sort_keys=True)

    print(f"\nSummary: Successfully registered {len(data_table_entries)} model(s)")


if __name__ == "__main__":
    main()
