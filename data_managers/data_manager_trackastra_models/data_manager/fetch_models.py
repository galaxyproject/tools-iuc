#!/usr/bin/env python3
"""Install versioned Trackastra model archives into a Galaxy data table."""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
import urllib.request
import zipfile
from hashlib import sha256
from pathlib import Path, PurePosixPath

DATA_TABLE_NAME = "trackastra_models"
REQUIRED_FILES = {"model.pt", "config.yaml", "train_config.yaml"}
MODEL_CATALOG = {
    "general_2d": {
        "name": "general_2d",
        "dimensionality": "2D",
        "version": "0.3.0",
        "url": (
            "https://github.com/weigertlab/trackastra-models/releases/"
            "download/v0.3.0/general_2d.zip"
        ),
    },
    "ctc": {
        "name": "ctc",
        "dimensionality": "2D/3D",
        "version": "0.3.0",
        "url": (
            "https://github.com/weigertlab/trackastra-models/releases/"
            "download/v0.3.0/ctc.zip"
        ),
    },
}


def safe_extract_zip(archive: Path, destination: Path) -> None:
    """Extract a ZIP without accepting absolute or parent-traversal paths."""
    destination = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            path = PurePosixPath(info.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Unsafe archive member: {info.filename}")
            target = (destination / Path(*path.parts)).resolve()
            if destination not in target.parents and target != destination:
                raise ValueError(f"Unsafe archive member: {info.filename}")
        handle.extractall(destination)


def find_model_directory(root: Path) -> Path:
    """Find the unique directory containing a complete Trackastra model."""
    candidates = [
        path
        for path in [root, *root.rglob("*")]
        if path.is_dir()
        and REQUIRED_FILES.issubset({item.name for item in path.iterdir()})
    ]
    if len(candidates) != 1:
        raise ValueError(
            "Expected exactly one directory containing model.pt, config.yaml, "
            f"and train_config.yaml; found {len(candidates)}."
        )
    return candidates[0]


def download(url: str, output: Path) -> None:
    """Download a model archive with a descriptive user agent."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Galaxy-Trackastra-data-manager"},
    )
    with (
        urllib.request.urlopen(request) as response,
        output.open("wb") as handle,
    ):
        shutil.copyfileobj(response, handle)


def install_archive(
    archive: Path,
    output_directory: Path,
    value: str,
) -> tuple[Path, str]:
    """Validate and install a model archive.

    Return the installed directory and archive hash.
    """
    archive_hash = sha256(archive.read_bytes()).hexdigest()
    destination = output_directory / value
    if destination.exists():
        raise ValueError(f"Model destination already exists: {destination}")

    with tempfile.TemporaryDirectory(prefix="trackastra-model-") as temp:
        extracted = Path(temp) / "extracted"
        extracted.mkdir()
        safe_extract_zip(archive, extracted)
        model_directory = find_model_directory(extracted)
        shutil.copytree(model_directory, destination)
    return destination, archive_hash


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_manager_json", type=Path)
    parser.add_argument("--known-models", default="")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--models",
        help="Comma-separated keys from the built-in catalog",
    )
    source.add_argument(
        "--archive",
        type=Path,
        help="Administrator-provided model ZIP",
    )
    parser.add_argument("--value")
    parser.add_argument("--dimensionality", choices=("2D", "3D", "2D/3D"))
    parser.add_argument("--model-version")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.data_manager_json.open() as handle:
        config = json.load(handle)
    try:
        output_directory = Path(config["output_data"][0]["extra_files_path"])
    except (KeyError, IndexError) as exc:
        raise ValueError(
            "Data manager output is missing extra_files_path."
        ) from exc
    output_directory.mkdir(parents=True, exist_ok=True)

    known_models = {item for item in args.known_models.split(",") if item}
    rows: list[dict[str, str]] = []

    if args.models:
        requested = [item for item in args.models.split(",") if item]
        for value in requested:
            if value not in MODEL_CATALOG:
                raise ValueError(f"Unknown Trackastra model: {value}")
            if value in known_models:
                raise ValueError(
                    f"Trackastra model is already installed: {value}"
                )
            metadata = MODEL_CATALOG[value]
            with tempfile.TemporaryDirectory(
                prefix="trackastra-download-"
            ) as temp:
                archive = Path(temp) / f"{value}.zip"
                download(metadata["url"], archive)
                installed, archive_hash = install_archive(
                    archive,
                    output_directory,
                    value,
                )
            rows.append(
                {
                    "value": value,
                    "name": metadata["name"],
                    "dimensionality": metadata["dimensionality"],
                    "version": metadata["version"],
                    "sha256": archive_hash,
                    "path": str(installed),
                    "source": metadata["url"],
                }
            )
    else:
        if not all((args.value, args.dimensionality, args.model_version)):
            raise ValueError(
                "Custom model installation requires value, dimensionality, "
                "and model version."
            )
        if args.value in known_models:
            raise ValueError(
                f"Trackastra model is already installed: {args.value}"
            )
        installed, archive_hash = install_archive(
            args.archive,
            output_directory,
            args.value,
        )
        rows.append(
            {
                "value": args.value,
                "name": args.value,
                "dimensionality": args.dimensionality,
                "version": args.model_version,
                "sha256": archive_hash,
                "path": str(installed),
                "source": "administrator_upload",
            }
        )

    result = {"data_tables": config.get("data_tables", {})}
    result["data_tables"][DATA_TABLE_NAME] = rows
    with args.data_manager_json.open("w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
