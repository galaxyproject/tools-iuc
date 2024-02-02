import argparse
import json
import os
import subprocess
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


def download_compleasm_database(name, output_dir):
    db_path = os.path.join(output_dir, f"{name}.db")
    if os.path.exists(db_path):
        print(f"Database {name} already exists in {output_dir}, skipping download.")
    else:
        command = f"compleasm download {name}"
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {name}: {e}")
            exit(1)
    return db_path


def main(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    data_manager_entries = []

    for name in args.name.split(','):
        if name not in busco_lineages:
            print(f"Error: Lineage '{name}' not found in available lineages.")
            exit(1)
        db_path = download_compleasm_database(name, args.output_dir)

        base_name = os.path.splitext(os.path.basename(db_path))[0]
        entry = {
            "value": base_name,
            "name": name,  # Utilisez le nom du lineage actuel ici
            "path": str(Path(args.output_dir)),
            "version": "0.2.5",
        }
        data_manager_entries.append(entry)

    data_manager_json = {"data_tables": {"compleasm_data": data_manager_entries}}

    with open(args.json, "w") as fh:
        json.dump(data_manager_json, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download data for CompleASM')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Output directory for saving databases')
    parser.add_argument('--name', required=True, help='Comma-separated list of lineage names to download')
    parser.add_argument('--json', required=True, help='Path to JSON file')
    args = parser.parse_args()
    main(args)
