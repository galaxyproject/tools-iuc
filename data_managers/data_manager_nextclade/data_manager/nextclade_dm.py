#!/usr/bin/env python

import argparse
import datetime
import itertools
import json
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


def get_dataset_list(
    names: List[str] = None,
    start_date: datetime.datetime = None,
    end_date: datetime.datetime = None,
    latest_only: bool = True
) -> List[dict]:
    """Retrieve and return a list of nextclade datasets

    When called with no arguments, returns the latest dataset of each dataset
    name/reference combination, equivalent to the result of running
    nextclade dataset list
    with no arguments.

    The list of datasets can be restricted by providing a list of dataset
    names.
    If an end_date is specified, returns the latest datasets published before
    or on this date.
    If latest_only is set to False, returns all datasets (filtered by the list
    of dataset names if specified) ever published or published within the
    specified time window.

    The return value is a list of dictionaries, which are themselves simplified
    versions (holding only selected key dataset info) of the json
    representation of the correpsonding datasets as obtained by calling
    nextclade dataset list --json
    on the commend line.
    """
    list_cmd = [
        "nextclade",
        "dataset",
        "list",
        "--json",
        "--include-old",
        "--include-incompatible",
    ]
    list_proc = subprocess.run(list_cmd, capture_output=True, check=True)
    dataset_list = sorted(
        json.loads(list_proc.stdout),
        key=lambda x: (
            x['attributes']['name']['value'],
            x['attributes']['reference']['value'],
            x['attributes']['tag']['value']),
        reverse=True
    )

    entry_list = []
    for dataset_key, versions in itertools.groupby(
        dataset_list,
        key=lambda x: (
            x['attributes']['name']['value'],
            x['attributes']['reference']['value']
        )
    ):
        name, reference = dataset_key
        if not names or name in names:
            for version in versions:
                tag = version["attributes"]["tag"]["value"]
                date = datetime.datetime.fromisoformat(tag.replace("Z", ""))
                if start_date and date < start_date:
                    break
                if not end_date or date <= end_date:
                    entry = {
                        "value": f"{name}_{reference}_{tag}",
                        "dataset_name": name,
                        "reference": reference,
                        "description": f"{tag}/{reference}",
                        "date": date,
                        "tag": tag,
                        "min_nextclade_version": version["compatibility"]["nextcladeCli"]["min"],
                        "max_nextclade_version": version["compatibility"]["nextcladeCli"]["max"] or ""
                    }
                    entry_list.append(entry)
                    if latest_only:
                        break
    return entry_list


def download(release: dict) -> pathlib.Path:
    """Download a nextclade dataset

    Given a simplified dictionary representation of a dataset as returned
    by get_dataset_list, downloads the dataset using the
    nextclade dataset get
    command.
    """
    dataset_name = release["dataset_name"]
    reference = release["reference"]
    tag = release["tag"]
    tag_file_suffix = tag.replace(":", "-")
    output_path = pathlib.Path(
        output_directory
    ) / f"{dataset}_{reference}_{tag_file_suffix}"

    download_cmd = [
        "nextclade",
        "dataset",
        "get",
        "--name",
        dataset_name,
        "--reference",
        reference,
        "--tag",
        tag,
        "--output-dir",
        output_path,
    ]
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
    parser.add_argument("datatable_name", default="nextclade2")
    parser.add_argument("galaxy_config")
    args = parser.parse_args()

    # known-revisions is populated from the Galaxy data table by the wrapper
    if args.known_revisions is not None:
        existing_release_tags = set(args.known_revisions)
    else:
        existing_release_tags = set()

    releases = [
        release for release in get_dataset_list(
            args.datasets, args.start_date, args.end_date, args.latest
        ) if release["value"] not in existing_release_tags
    ]

    with open(args.galaxy_config) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)

    data_manager_dict = {"data_tables": {args.datatable_name: []}}

    for release in releases:
        if args.testmode:
            for release in releases:
                for k, v in release.items():
                    print(k, v)
        else:
            fname = download(release)
            if fname is not None:
                data_manager_dict["data_tables"][args.datatable_name].append(
                    {
                        "value": release["value"],
                        "dataset_name": release["dataset_name"],
                        "description": release["description"],
                        # date in ISO 8601 is a good key to sort on
                        "date": release["date"].isoformat(),
                        # record compatible nextclade versions for filtering
                        # max_nextclade_version will typically be an empty
                        # string for recent datasets at install time, but an
                        # admin could update the value later to prevent the
                        # record from showing up in the interface of more
                        # recent, incompatible nextclade wrapper versions.
                        "min_nextclade_version": release["min_nextclade_version"],
                        "max_nextclade_version": release["max_nextclade_version"],
                        "path": str(output_directory / fname),
                    }
                )

    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=lambda x: x["value"], reverse=True
    )
    with open(args.galaxy_config, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
