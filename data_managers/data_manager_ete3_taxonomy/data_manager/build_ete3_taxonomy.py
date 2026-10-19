import argparse
import hashlib
import tarfile
import tempfile
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


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output-db",
        required=True,
        help="Path where the ETE3 SQLite database will be written.",
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
        "--identity-output",
        required=True,
        help="File in which to write the stable taxonomy SHA-256 digest.",
    )

    args = parser.parse_args()

    output_db = Path(args.output_db)
    identity_output = Path(args.identity_output)

    output_db.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.taxonomy_dir:
        source_digest = build_from_existing_taxonomy(
            output_db,
            args.taxonomy_dir,
        )
    else:
        source_digest = build_from_download(
            output_db,
        )

    if not output_db.exists():
        raise RuntimeError(
            f"ETE3 taxonomy database was not created: {output_db}"
        )

    if output_db.stat().st_size == 0:
        raise RuntimeError(
            f"ETE3 taxonomy database is empty: {output_db}"
        )

    identity_output.write_text(
        f"{source_digest}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
