import argparse
import tarfile
import tempfile
from pathlib import Path

from ete3 import NCBITaxa
from ete3.ncbi_taxonomy import ncbiquery


def build_from_existing_taxonomy(output_db, taxonomy_dir):
    taxonomy_dir = Path(taxonomy_dir)

    required = [
        "nodes.dmp",
        "names.dmp",
    ]

    for filename in required:
        path = taxonomy_dir / filename
        if not path.exists():
            raise RuntimeError(
                f"Required NCBI taxonomy file not found: {path}"
            )

    with tempfile.NamedTemporaryFile(
        suffix=".tar.gz",
        delete=False,
    ) as handle:
        taxdump = Path(handle.name)

    with tarfile.open(taxdump, "w:gz") as archive:
        for path in taxonomy_dir.iterdir():
            if path.is_file():
                archive.add(path, arcname=path.name)

    NCBITaxa(
        dbfile=str(output_db),
        taxdump_file=str(taxdump),
    )


def build_from_download(output_db):
    # ETE3's update_db downloads the NCBI taxdump when no
    # taxdump_file is supplied.
    ncbiquery.update_db(str(output_db))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output-db",
        required=True,
    )

    parser.add_argument(
        "--taxonomy-dir",
        default=None,
    )

    args = parser.parse_args()

    output_db = Path(args.output_db)

    if args.taxonomy_dir:
        build_from_existing_taxonomy(
            output_db,
            args.taxonomy_dir,
        )
    else:
        build_from_download(output_db)

    if not output_db.exists():
        raise RuntimeError(
            f"ETE3 taxonomy database was not created: {output_db}"
        )


if __name__ == "__main__":
    main()