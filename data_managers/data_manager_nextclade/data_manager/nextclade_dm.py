#!/usr/bin/env python

import argparse
import datetime
import json
import operator
import pathlib
import subprocess
import sys
from typing import List


def parse_date(d: str) -> datetime.datetime:
    # Parses the publication date from the nextclade release tags or user input into a datetime object.
    date = None
    try:
        date = datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        date = datetime.datetime.strptime(d, "%Y-%m-%d")
    return date


def entry_to_tag(entry: dict) -> str:
    return (
        entry["attributes"]["name"]["value"] + "_" + entry["attributes"]["tag"]["value"]
    )


def get_database_list() -> List[dict]:
    list_cmd = [
        "nextclade",
        "dataset",
        "list",
        "--json",
        "--include-old",
        "--include-incompatible",
    ]
    list_proc = subprocess.run(list_cmd, capture_output=True, check=True)
    database_list = json.loads(list_proc.stdout)
    entry_list = []
    for db_entry in database_list:
        attributes = db_entry["attributes"]
        entry = {
            "value": entry_to_tag(db_entry),
            "database_name": attributes["name"]["value"],
            "description": attributes["name"]["valueFriendly"],
            "date": datetime.datetime.fromisoformat(
                attributes["tag"]["value"].replace("Z", "")
            ),
            "tag": attributes["tag"]["value"],
            "min_nextclade_version": db_entry["compatibility"]["nextcladeCli"]["min"],
        }
        entry_list.append(entry)
    return entry_list


def filter_by_date(
    existing_release_tags: List[str],
    name: str,
    releases: list,
    start_date: datetime.datetime = None,
    end_date: datetime.datetime = None,
) -> List[dict]:
    ret = []
    for release in releases:
        if (
            release["database_name"] != name
            or release["value"] in existing_release_tags
        ):
            continue
        if start_date and release["date"] < start_date:
            break
        if not end_date or release["date"] <= end_date:
            ret.append(release)

    return ret


def download_and_unpack(name: str, release: str, output_directory: str) -> pathlib.Path:
    download_cmd = [
        "nextclade",
        "dataset",
        "get",
        "--name",
        name,
        "--tag",
        release,
        "--output-dir",
    ]
    output_path = pathlib.Path(output_directory) / (
        name + "_" + release.replace(":", "-")
    )
    download_cmd.append(str(output_path))
    subprocess.run(download_cmd, check=True)
    return output_path


def comma_split(args: str) -> List[str]:
    return args.split(",")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--testmode", default=False, action="store_true")
    parser.add_argument("--latest", default=False, action="store_true")
    parser.add_argument("--start_date", type=parse_date)
    parser.add_argument("--end_date", type=parse_date)
    parser.add_argument("--known_revisions", type=comma_split)
    parser.add_argument("--datasets", type=comma_split, default=["sars-cov-2"])
    parser.add_argument("datatable_name", default="nextclade")
    parser.add_argument("galaxy_config")
    args = parser.parse_args()

    # known-revisions is populated from the Galaxy data table by the wrapper
    if args.known_revisions is not None:
        existing_release_tags = set(args.known_revisions)
    else:
        existing_release_tags = set()

    releases_available = get_database_list()
    if args.testmode:
        releases = []
        for name in args.datasets:
            releases.extend(
                filter_by_date(
                    [],
                    name,
                    releases_available,
                    start_date=args.start_date,
                    end_date=args.end_date,
                )
            )
        for release in releases:
            print(
                release["value"],
                release["description"],
                release["date"].isoformat(),
                release["min_nextclade_version"],
            )
        sys.exit(0)

    with open(args.galaxy_config) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)

    data_manager_dict = {"data_tables": {args.datatable_name: []}}

    releases = []
    if args.latest:
        for dataset in args.datasets:
            for release in releases_available:
                if release["database_name"] == dataset:
                    if release["value"] not in existing_release_tags:
                        # add the latest release for this dataset, but only if we don't already have it
                        releases.append(release)
                    break
    else:
        for dataset in args.datasets:
            releases_for_ds = filter_by_date(
                existing_release_tags,
                dataset,
                releases_available,
                start_date=args.start_date,
                end_date=args.end_date,
            )
            releases.extend(releases_for_ds)

    for release in releases:
        fname = download_and_unpack(
            release["database_name"], release["tag"], output_directory
        )
        if fname is not None:
            data_manager_dict["data_tables"][args.datatable_name].append(
                {
                    "value": release["value"],
                    "database_name": release["database_name"],
                    "description": release["description"],
                    "min_nextclade_version": release["min_nextclade_version"],
                    "date": release["date"].isoformat(),  # ISO 8601 is easily sortable
                    "path": str(output_directory / fname),
                }
            )
    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_config, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
