import argparse
import json


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--output", required=True)
    parser.add_argument("--path", required=True)
    parser.add_argument("--value", required=True)
    parser.add_argument("--name", required=True)

    args = parser.parse_args()

    data_manager_json = {
        "data_tables": {
            "ete3_taxonomy_db": [
                {
                    "value": args.value,
                    "name": args.name,
                    "path": args.path,
                }
            ]
        }
    }

    with open(args.output, "w") as handle:
        json.dump(data_manager_json, handle)


if __name__ == "__main__":
    main()