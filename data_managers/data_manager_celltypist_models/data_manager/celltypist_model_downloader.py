import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

MODEL_JSON_URL = "https://celltypist.cog.sanger.ac.uk/models/models.json"


def fetch_json(url):
    with urllib.request.urlopen(url) as r:
        return json.load(r)


def safe_download(url, dest):
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} to {dest}")
    urllib.request.urlretrieve(url, dest)


def model_version(filename, version):
    base = Path(filename).stem.replace(" ", "_").replace("/", "_")
    return f"{base}_{version}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, help="Output JSON file")
    p.add_argument("--model", default="all", help="Model name to download, or 'all' for all models")
    p.add_argument("--version", default="", help="Model version to download, or empty for latest")
    args = p.parse_args()

    with open(args.out) as fh:
        params = json.load(fh)
    target_directory = Path(params["output_data"][0]["extra_files_path"])
    target_directory.mkdir(parents=True, exist_ok=True)

    print(f"Fetching model metadata from {MODEL_JSON_URL}...", file=sys.stderr)
    models = fetch_json(MODEL_JSON_URL).get("models", [])

    # Filter by model name if specified
    if args.model and args.model != "all":
        target = args.model
        models = [m for m in models if m.get("filename", "").replace(".pkl", "") == target]

    # Handle version filtering
    if args.version and models:
        requested_version = args.version if args.version.startswith("v") else f"v{args.version}"

        versioned_model = next((m for m in models if m.get("version", "") == requested_version), None)
        if versioned_model:
            models = [versioned_model]
        else:
            latest_model = models[0]
            latest_url = latest_model.get("url", "")
            if latest_url:
                new_url = re.sub(r"/v\d+/", f"/{requested_version}/", latest_url)
                if new_url != latest_url:
                    modified_model = latest_model.copy()
                    modified_model["url"] = new_url
                    modified_model["version"] = requested_version
                    models = [modified_model]
                    print(
                        f"Warning: Version {requested_version} not in catalog. Attempting to use URL: {new_url}",
                        file=sys.stderr,
                    )
    # If no version specified, use all models from JSON (which contains only latest versions)

    out_json = {"data_tables": {"celltypist_models": []}}

    for m in models:
        filename = m.get("filename", "")
        url = m.get("url")
        name = m.get("name", "")
        version = m.get("version", "")
        date = m.get("date", "")

        if name and version:
            name_with_version = f"{name} ({version})"
        else:
            name_with_version = name

        model = model_version(filename, version)

        # Download model files
        if url:
            dest_path = target_directory / f"{model}.pkl"
            if not dest_path.exists():
                print(f"Downloading {filename}...", file=sys.stderr)
                safe_download(url, str(dest_path))
            path_value = str(dest_path)
        else:
            path_value = ""

        out_json["data_tables"]["celltypist_models"].append({
            "value": model,
            "name": name_with_version,
            "date": date,
            "path": path_value,
        })

    # Write back the data_manager_json with updated data tables
    with open(args.out, "w") as f:
        json.dump(out_json, f, indent=2)


if __name__ == "__main__":
    main()
