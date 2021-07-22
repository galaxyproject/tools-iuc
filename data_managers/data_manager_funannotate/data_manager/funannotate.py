#!/usr/bin/env python

import argparse
import json
import operator
import os
import subprocess
import sys

from datetime import datetime


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--partial', dest='partial', action='store_true', help='Only download a small subset of data (for testing)')
    parser.add_argument("version_id")
    parser.add_argument("datatable_name")
    parser.add_argument("galaxy_datamanager_filename")
    args = parser.parse_args()

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][args.datatable_name] = data_manager_dict[
        "data_tables"
    ].get(args.datatable_name, [])

    os.mkdir(output_directory)
    cmd_args = ['funannotate', 'setup', '-d', output_directory, '-b', 'all']
    if args.partial:
        cmd_args += ['-i', 'merops', '-b', 'eukaryota']
    proc = subprocess.Popen(args=cmd_args, shell=False, cwd=output_directory)
    return_code = proc.wait()
    if return_code:
        print("Error downloading Funannotate database.", file=sys.stderr)
        sys.exit(return_code)

    version_id = datetime.today().strftime('%Y-%m-%d-%H%M%S')

    version = '1.0'

    data_manager_dict["data_tables"][args.datatable_name].append(
        dict(
            value=version_id,
            description="Funannotate database %s" % version_id,
            format_version=version,
            path=output_directory,
        )
    )

    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
