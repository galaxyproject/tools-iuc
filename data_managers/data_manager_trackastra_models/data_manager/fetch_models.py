#!/usr/bin/env python3
"""Install Trackastra model archives into a Galaxy tool-data table."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import stat
import tempfile
import urllib.request
import zipfile
from hashlib import sha256
from pathlib import Path, PurePosixPath
from typing import BinaryIO

DATA_TABLE_NAME = "trackastra_models"
IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")
REQUIRED_FILES = {"model.pt", "config.yaml", "train_config.yaml"}
DOWNLOAD_TIMEOUT = 120
MODEL_CATALOG = {
    "general_2d_0.3.0": {
        "name": "general_2d (0.3.0)",
        "dimensionality": "2D",
        "version": "0.3.0",
        "url": (
            "https://github.com/weigertlab/trackastra-models/releases/"
            "download/v0.3.0/general_2d.zip"
        ),
    },
    "ctc_0.3.0": {
        "name": "ctc (0.3.0)",
        "dimensionality": "2D/3D",
        "version": "0.3.0",
        "url": (
            "https://github.com/weigertlab/trackastra-models/releases/"
            "download/v0.3.0/ctc.zip"
        ),
    },
}


def validate_identifier(value: str, field_name: str) -> None:
    if not IDENTIFIER_PATTERN.fullmatch(value):
        raise ValueError(
            f"Invalid {field_name} {value!r}. Use only letters, numbers, dot, "
            "underscore, plus, and hyphen, beginning with a letter or number."
        )


def stream_sha256(handle: BinaryIO) -> str:
    digest = sha256()
    while chunk := handle.read(1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()


def archive_sha256(archive: Path) -> str:
    with archive.open("rb") as handle:
        return stream_sha256(handle)


def safe_extract_zip(archive: Path, destination: Path) -> None:
    """Extract a ZIP while rejecting traversal, links, and encrypted members."""
    if not zipfile.is_zipfile(archive):
        raise ValueError(f"Not a valid ZIP archive: {archive}")

    destination = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            path = PurePosixPath(info.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Unsafe archive member: {info.filename}")
            if info.flag_bits & 0x1:
                raise ValueError(
                    f"Encrypted archive members are not supported: {info.filename}"
                )
            member_mode = info.external_attr >> 16
            if stat.S_ISLNK(member_mode):
                raise ValueError(f"Symbolic links are not allowed: {info.filename}")

            target = (destination / Path(*path.parts)).resolve()
            if target != destination and destination not in target.parents:
                raise ValueError(f"Unsafe archive member: {info.filename}")
        handle.extractall(destination)


def find_model_directory(root: Path) -> Path:
    """Return the unique directory containing a complete Trackastra model."""
    candidates: list[Path] = []
    for path in (root, *root.rglob("*")):
        if not path.is_dir():
            continue
        children = {item.name: item for item in path.iterdir()}
        if not REQUIRED_FILES.issubset(children):
            continue
        if not all(children[name].is_file() for name in REQUIRED_FILES):
            raise ValueError(
                f"Required model entries must be regular files in {path}."
            )
        candidates.append(path)

    if len(candidates) != 1:
        raise ValueError(
            "Expected exactly one directory containing model.pt, config.yaml, "
            f"and train_config.yaml; found {len(candidates)}."
        )
    return candidates[0]


def download(url: str, output: Path) -> None:
    """Download one model archive."""
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Galaxy-Trackastra-data-manager"},
    )
    with (
        urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT) as response,
        output.open("wb") as handle,
    ):
        shutil.copyfileobj(response, handle)
    if output.stat().st_size == 0:
        raise ValueError(f"Downloaded an empty model archive from {url}")


def install_archive(
    archive: Path,
    output_directory: Path,
    value: str,
) -> tuple[Path, str]:
    """Validate and copy one model archive into the data-manager output."""
    validate_identifier(value, "model identifier")
    archive_hash = archive_sha256(archive)
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
        "--model",
        action="append",
        dest="models",
        help="Identifier from the built-in model catalog; may be repeated",
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


def model_row(
    *,
    value: str,
    name: str,
    dimensionality: str,
    version: str,
    archive_hash: str,
    installed: Path,
    source: str,
) -> dict[str, str]:
    return {
        "value": value,
        "name": name,
        "dimensionality": dimensionality,
        "version": version,
        "sha256": archive_hash,
        "path": str(installed),
        "source": source,
    }


def main() -> None:
    args = parse_args()
    with args.data_manager_json.open() as handle:
        config = json.load(handle)

    try:
        output_directory = Path(config["output_data"][0]["extra_files_path"])
    except (KeyError, IndexError, TypeError) as exc:
        raise ValueError(
            "Data manager output is missing extra_files_path."
        ) from exc
    output_directory.mkdir(parents=True, exist_ok=True)

    known_models = {item for item in args.known_models.split(",") if item}
    rows: list[dict[str, str]] = []

    if args.models is not None:
        requested = [item for item in args.models if item]
        if not requested:
            raise ValueError("Select at least one Trackastra model to install.")
        if len(requested) != len(set(requested)):
            raise ValueError("The model selection contains duplicate identifiers.")

        for value in requested:
            if value not in MODEL_CATALOG:
                raise ValueError(f"Unknown Trackastra model: {value}")
            if value in known_models:
                raise ValueError(f"Trackastra model is already installed: {value}")

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
                model_row(
                    value=value,
                    name=metadata["name"],
                    dimensionality=metadata["dimensionality"],
                    version=metadata["version"],
                    archive_hash=archive_hash,
                    installed=installed,
                    source=metadata["url"],
                )
            )
    else:
        if args.archive is None:
            raise ValueError("A custom model archive is required.")
        if not all((args.value, args.dimensionality, args.model_version)):
            raise ValueError(
                "Custom model installation requires value, dimensionality, "
                "and model version."
            )
        validate_identifier(args.value, "model identifier")
        validate_identifier(args.model_version, "model version")
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
            model_row(
                value=args.value,
                name=args.value,
                dimensionality=args.dimensionality,
                version=args.model_version,
                archive_hash=archive_hash,
                installed=installed,
                source="administrator_upload",
            )
        )

    result = {"data_tables": config.get("data_tables", {})}
    result["data_tables"][DATA_TABLE_NAME] = rows
    with args.data_manager_json.open("w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
