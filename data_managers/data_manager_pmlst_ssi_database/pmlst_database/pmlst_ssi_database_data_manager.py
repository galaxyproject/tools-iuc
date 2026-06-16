#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

PMLST_TABLE = "pmlst_db"


def database_key(version: str) -> str:
    safe_version = "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in version
    ).strip("._-")
    return f"pmlst_{safe_version or 'unknown'}"


def version_from_database(database: Path) -> str:
    version_file = database / "version.txt"
    if not version_file.exists():
        return "unknown"

    for line in version_file.read_text(encoding="utf-8").splitlines():
        if line.startswith("version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def install_database(destination: Path, fixture_db: Path | None) -> None:
    if fixture_db is not None:
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(fixture_db, destination)
        return

    subprocess.run(
        ["pmlst-download-db", str(destination)],
        check=True,
    )


def data_manager_row(database: Path) -> dict[str, str]:
    version = version_from_database(database)
    return {
        "value": database_key(version),
        "name": version,
        "path": str(database),
        "version": version,
    }


def write_data_manager_json(
    out_file: Path,
    row: dict[str, str] | None,
) -> None:
    data_manager_json: dict[str, Any] = {
        "data_tables": {
            PMLST_TABLE: [row] if row is not None else [],
        }
    }
    out_file.write_text(json.dumps(data_manager_json, sort_keys=True), encoding="utf-8")


def managed_database_target(
    data_manager_path: Path,
    value: str,
) -> Path:
    base_path = data_manager_path.resolve()
    target = (base_path / PMLST_TABLE / value).resolve()
    if not target.is_relative_to(base_path):
        raise ValueError(
            f"Refusing to remove target outside data-manager path: {target}"
        )
    return target


def remove_temporary_database(destination: Path, managed_target: Path) -> None:
    if destination.resolve() == managed_target.resolve():
        return
    if destination.is_dir():
        shutil.rmtree(destination)
    elif destination.exists():
        destination.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Install a pMLST database for Galaxy.")
    parser.add_argument("--out-file", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--data-manager-path", required=True, type=Path)
    parser.add_argument(
        "--fixture-db",
        type=Path,
        help="Copy this local fixture database instead of downloading upstream data.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    destination = args.destination.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fixture_db = args.fixture_db.resolve() if args.fixture_db else None
    install_database(destination, fixture_db)
    row = data_manager_row(destination)
    value = row["value"]
    managed_target = managed_database_target(args.data_manager_path.resolve(), value)
    if managed_target.exists():
        remove_temporary_database(destination, managed_target)
        write_data_manager_json(args.out_file, None)
    else:
        write_data_manager_json(args.out_file, row)


if __name__ == "__main__":
    main()
