import argparse
import hashlib
import json
import sqlite3
import tarfile
import tempfile
from datetime import date
from pathlib import Path
from urllib.request import urlretrieve

from ete3 import NCBITaxa


NCBI_TAXDUMP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
)


def taxonomy_digest(taxonomy_dir):
    """Return a deterministic SHA-256 digest of a taxonomy directory."""
    taxonomy_dir = Path(taxonomy_dir)

    manifest = []

    for path in sorted(taxonomy_dir.iterdir()):
        if not path.is_file():
            continue

        file_hash = hashlib.sha256()

        with path.open("rb") as handle:
            for chunk in iter(
                lambda: handle.read(1024 * 1024),
                b"",
            ):
                file_hash.update(chunk)

        manifest.append(
            f"{path.name}\t{file_hash.hexdigest()}"
        )

    digest = hashlib.sha256()

    for entry in manifest:
        digest.update(entry.encode("utf-8"))
        digest.update(b"\n")

    return digest.hexdigest()


def validate_taxonomy_dir(taxonomy_dir):
    """Check that the files required by ETE3 are present."""
    taxonomy_dir = Path(taxonomy_dir)

    required = [
        "nodes.dmp",
        "names.dmp",
        "merged.dmp",
    ]

    for filename in required:
        path = taxonomy_dir / filename

        if not path.exists():
            raise RuntimeError(
                f"Required NCBI taxonomy file not found: {path}"
            )


def build_from_existing_taxonomy(output_db, taxonomy_dir):
    """Build the ETE3 database from an installed NCBI taxonomy."""
    taxonomy_dir = Path(taxonomy_dir)

    validate_taxonomy_dir(taxonomy_dir)

    source_digest = taxonomy_digest(taxonomy_dir)

    with tempfile.TemporaryDirectory() as tmpdir:
        taxdump = Path(tmpdir) / "taxdump.tar.gz"

        with tarfile.open(taxdump, "w:gz") as archive:
            for path in sorted(taxonomy_dir.iterdir()):
                if path.is_file():
                    archive.add(
                        path,
                        arcname=path.name,
                    )

        NCBITaxa(
            dbfile=str(output_db),
            taxdump_file=str(taxdump),
        )

    return source_digest


def build_from_download(output_db):
    """Download the current NCBI taxdump and build the ETE3 database."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        taxdump = tmpdir / "taxdump.tar.gz"
        taxonomy_dir = tmpdir / "taxonomy"
        taxonomy_dir.mkdir()

        urlretrieve(
            NCBI_TAXDUMP_URL,
            taxdump,
        )

        with tarfile.open(taxdump, "r:gz") as archive:
            archive.extractall(
                taxonomy_dir,
                filter="data",
            )

        validate_taxonomy_dir(taxonomy_dir)

        source_digest = taxonomy_digest(taxonomy_dir)

        NCBITaxa(
            dbfile=str(output_db),
            taxdump_file=str(taxdump),
        )

    return source_digest


def database_version(output_db):
    """Return the schema version ETE3 recorded in the database it built."""
    with sqlite3.connect(output_db) as db:
        row = db.execute("SELECT version FROM stats").fetchone()

    if row is None:
        raise RuntimeError(
            f"ETE3 taxonomy database records no schema version: {output_db}"
        )

    return str(row[0])


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output-dir",
        required=True,
        help=(
            "Directory (a direct child of the output's extra files path) that "
            "will hold taxa.sqlite. ETE3 also writes its taxa.sqlite.traverse.pkl "
            "cache next to the database, so the whole directory is registered in "
            "the data table."
        ),
    )

    parser.add_argument(
        "--output-json",
        required=True,
        help="Path where the JSON output will be written.",
    )

    parser.add_argument(
        "--taxonomy-dir",
        default=None,
        help=(
            "Existing directory containing the NCBI taxonomy dump. "
            "If omitted, the current NCBI taxdump is downloaded."
        ),
    )

    parser.add_argument(
        "--taxonomy-value",
        default=None,
        help="Value of the existing NCBI taxonomy dump.",
    )

    parser.add_argument(
        "--taxonomy-name",
        default=None,
        help="Name of the existing NCBI taxonomy dump.",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_db = output_dir / "taxa.sqlite"

    if args.taxonomy_dir:
        source_digest = build_from_existing_taxonomy(
            output_db,
            args.taxonomy_dir,
        )

        value = args.taxonomy_value
        name = args.taxonomy_name

    else:
        source_digest = build_from_download(output_db)

        value = f"ncbi_{source_digest}"
        name = f"NCBI Taxonomy (downloaded {date.today().isoformat()})"

    if not output_db.exists():
        raise RuntimeError(
            f"ETE3 taxonomy database was not created: {output_db}"
        )

    if output_db.stat().st_size == 0:
        raise RuntimeError(
            f"ETE3 taxonomy database is empty: {output_db}"
        )

    version = database_version(output_db)

    data_manager_json = {
        "data_tables": {
            "ete3_taxonomy_db": [
                {
                    "value": value,
                    "name": name,
                    "version": version,
                    "path": output_dir.name,
                }
            ]
        }
    }

    with open(args.output_json, "w") as handle:
        json.dump(data_manager_json, handle)


if __name__ == "__main__":
    main()
