import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--path", required=True)
    args = parser.parse_args()

    data_manager_json = {
        "data_tables": {
            "perseus_taxonomy_db": [
                {
                    "value": "ncbi_taxonomy",
                    "name": "NCBI Taxonomy",
                    "path": args.path,
                }
            ]
        }
    }

    with open(args.output, "w") as handle:
        json.dump(data_manager_json, handle)


if __name__ == "__main__":
    main()
