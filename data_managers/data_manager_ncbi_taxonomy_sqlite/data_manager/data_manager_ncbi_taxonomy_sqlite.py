from __future__ import division, print_function

import argparse
import datetime
import json
import os
import os.path
import shlex
import subprocess

DATA_TABLE_NAME = "ncbi_taxonomy_sqlite"


def build_sqlite(taxonomy_dir, output_directory, name=None, description=None):
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    output_filename = os.path.join(output_directory, "tax.ncbitaxonomy.sqlite")
    cmd_str = "taxonomy_util -d {} to_sqlite {}".format(output_filename, taxonomy_dir)
    cmd = shlex.split(cmd_str)
    subprocess.check_call(cmd)

    today_str = datetime.date.today().strftime("%Y-%m-%d")
    if name is None or name.strip() == "":
        name = "ncbitaxonomy_build_" + today_str

    if description is None or description.strip() == "":
        description = "NCBI Taxonomy database (built on {})".format(today_str)

    data = [dict(value=name, description=description, path=output_filename)]
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build SQLite database from NCBI taxonomy"
    )
    parser.add_argument(
        "--output_directory", default="tmp", help="Directory to write output to"
    )
    parser.add_argument(
        "taxonomy_dir",
        help="Path to directory containing NCBI Taxonomy nodes.dml and names.dmp file"
    )
    parser.add_argument(
        "name",
        help="Name to use for the entry in the data table"
    )
    parser.add_argument(
        "description",
        help="Description to use for the entry in the data table"
    )
    parser.add_argument(
        "galaxy_datamanager_filename",
        help="Galaxy JSON format file describing data manager inputs",
    )
    args = parser.parse_args()

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)
    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    if output_directory is None:
        output_directory = args.output_directory

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = data_manager_dict[
        "data_tables"
    ].get(DATA_TABLE_NAME, [])

    data = build_sqlite(args.taxonomy_dir, output_directory, args.name, args.description)

    data_manager_dict["data_tables"][DATA_TABLE_NAME].extend(data)
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)
