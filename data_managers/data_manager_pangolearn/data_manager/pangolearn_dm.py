#!/usr/bin/env python

import argparse
import datetime
import json
import operator
import os
import shutil
import sys
import tarfile

import requests


def extract_date(tag_str):
    parts = tag_str.split("_")
    assert len(parts) < 3, "expected maximum of two parts, got " + str(parts)
    tag_date = datetime.datetime.strptime(parts[0], "%Y-%m-%d")
    if len(parts) == 2:
        version = int(parts[1])
        assert (
            version < 24 * 60
        )  # because the code stores versions as minutes of the day, it can't handle versions > 1440
        tag_date += datetime.timedelta(minutes=version)
    return tag_date


def get_model_list(
    existing_release_tags,
    url="https://api.github.com/repos/cov-lineages/pangoLEARN/releases",
):
    response = requests.get(url)
    if response.status_code == 200:
        release_list = json.loads(response.text)
        release_info = [
            dict(
                tag_name=e["tag_name"],
                name=e["name"],
                date=extract_date(e["tag_name"]),
                tarball_url=e["tarball_url"],
            )
            for e in release_list
            if e["tag_name"] not in existing_release_tags
        ]
        return release_info
    else:
        response.raise_for_status()


def filter_by_date(existing_release_tags, start_date=None, end_date=None):
    release_list = get_model_list(existing_release_tags)
    return [
        element
        for element in release_list
        if not (
            (end_date is not None and element["date"] > end_date)
            or (start_date is not None and element["date"] < start_date)
        )
    ]


def download_and_unpack(url, output_directory):
    response = requests.get(url)
    if response.status_code == 200:
        tmp_filename = url.split("/")[-1]
        tmpfile = open(tmp_filename, "wb")
        tmpfile.write(response.content)
        tmpfile.close()
        shutil.copy(tmp_filename, "/tmp")
        tf = tarfile.open(tmp_filename)
        pl_path = tf.next().name
        tf.extractall(output_directory)
        os.unlink(tmp_filename)
        os.rename(
            output_directory + "/" + pl_path + "/" + "pangoLEARN",
            output_directory + "/" + tmp_filename,
        )
        shutil.rmtree(output_directory + "/" + pl_path)
        return tmp_filename
    else:
        response.raise_for_status()


def parse_date(d):
    return datetime.datetime.strptime(d, "%Y-%m-%d")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--testmode", default=False, action="store_true")
    parser.add_argument("--latest", default=False, action="store_true")
    parser.add_argument("--start_date", type=parse_date)
    parser.add_argument("--end_date", type=parse_date)
    parser.add_argument("--overwrite", default=False, action="store_true")
    parser.add_argument('--pangolearn_format_version')
    parser.add_argument("datatable_name")
    parser.add_argument("galaxy_datamanager_filename")
    args = parser.parse_args()

    if args.testmode:
        releases = filter_by_date(start_date=args.start_date, end_date=args.end_date)
        for release in releases:
            print(release["tag_name"], release["tarball_url"].split("/")[-1])
        sys.exit(0)

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][args.datatable_name] = data_manager_dict[
        "data_tables"
    ].get(args.datatable_name, [])

    # NOTE: the data_manager_dict["data_tables"][args.datatable_name] is not actually populated with the
    # contents of the existing data table, so the "no-overwrite" logic and the
    # only-download-what-we-don't-have logic does not in fact work. It is left but unused for now.
    if not args.overwrite:
        existing_release_tags = set(
            [
                el["value"]
                for el in data_manager_dict["data_tables"][args.datatable_name]
            ]
        )
    else:
        existing_release_tags = set()
    if args.latest:
        releases = [get_model_list(existing_release_tags)[0]]
    else:
        releases = filter_by_date(
            existing_release_tags, start_date=args.start_date, end_date=args.end_date
        )
    releases_to_download = [
        release
        for release in releases
        if release["tag_name"] not in existing_release_tags
    ]
    for release in releases_to_download:
        tag = download_and_unpack(release["tarball_url"], output_directory)
        release_date = parse_date(tag)
        if args.pangolearn_format_version is not None:
            version = args.pangolearn_format_version
        else:
            # 2021-05-27 was the first release of pangoLEARN for pangolin 3, which changed DB format
            if release_date >= datetime.datetime(2021, 5, 27):
                version = '2.0'
            else:
                version = '1.0'
        data_manager_dict["data_tables"][args.datatable_name].append(
            dict(
                value=tag,
                description=release["name"],
                format_version=version,
                path=output_directory + "/" + tag,
            )
        )
    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
