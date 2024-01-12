import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

busco_lineages = [
    "aconoidasida", "actinopterygii", "agaricales", "agaricomycetes", "alveolata",
    "apicomplexa", "arachnida", "arthropoda", "ascomycota", "aves",
    "basidiomycota", "boletales", "brassicales", "capnodiales", "carnivora",
    "cetartiodactyla", "chaetothyriales", "chlorophyta", "coccidia", "cyprinodontiformes",
    "diptera", "dothideomycetes", "embryophyta", "endopterygota", "euarchontoglires",
    "eudicots", "euglenozoa", "eukaryota", "eurotiales",
    "eurotiomycetes", "eutheria", "fabales", "fungi", "glires",
    "glomerellales", "helotiales", "hemiptera", "hymenoptera", "hypocreales",
    "insecta", "laurasiatheria", "leotiomycetes", "lepidoptera", "liliopsida",
    "mammalia", "metazoa", "microsporidia", "mollusca", "mucorales",
    "mucoromycota", "nematoda", "onygenales", "passeriformes", "plasmodium",
    "pleosporales", "poales", "polyporales", "primates", "saccharomycetes",
    "sauropsida", "solanales", "sordariomycetes", "stramenopiles", "tetrapoda",
    "tremellomycetes", "vertebrata", "viridiplantae",
]


def download_compleasm_database(lineage, output_dir):
    command = f"compleasm download {lineage}"
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {lineage}: {e}")
        sys.exit(1)

    return os.path.join("mb_downloads", f"{lineage}.db")


def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    downloaded_databases = []
    for lineage in busco_lineages:
        db_path = download_compleasm_database(lineage, args.output_dir)
        downloaded_databases.append(db_path)

    data_manager_entries = []
    for db_path in downloaded_databases:
        base_name = os.path.splitext(os.path.basename(db_path))[0]
        entry = {
            "value": base_name,
            "name": args.name,
            "path": str(Path(args.output_dir)),
        }
        data_manager_entries.append(entry)

    data_manager_json = {"data_tables": {"compleasm_data": data_manager_entries}}

    with open(args.json, "w") as fh:
        json.dump(data_manager_json, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download data for CompleASM')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Output directory for saving databases')
    parser.add_argument('--name', default='default_name', help='Data table entry unique ID')
    parser.add_argument('--json', help='Path to JSON file')
    args = parser.parse_args()

    main(args)
