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


def get_model_list(
    existing_release_tags,
    url="https://api.github.com/repos/cov-lineages/pangoLEARN/releases"
):
    page_num = 0
    while True:
        page_num += 1
        response = requests.get(url + f'?page={page_num}')
        if response.status_code == 200:
            release_list_chunk = json.loads(response.text)
            if not release_list_chunk:
                # past the last page of results
                return
            for e in release_list_chunk:
                if e["tag_name"] in existing_release_tags:
                    continue
                if e["prerelease"]:
                    continue
                yield dict(
                    tag_name=e["tag_name"],
                    name=e["name"],
                    date=parse_date(e["tag_name"]),
                    tarball_url=e["tarball_url"],
                )
        else:
            response.raise_for_status()


def filter_by_date(existing_release_tags, start_date=None, end_date=None):
    ret = []
    for release in get_model_list(existing_release_tags):
        if start_date and release["date"] < start_date:
            break
        if not end_date or release["date"] <= end_date:
            ret.append(release)

    return ret


def download_and_unpack(url, output_directory, v3datatree=True):
    response = requests.get(url)
    if response.status_code == 200:
        tmp_filename = url.split("/")[-1]
        tmpfile = open(tmp_filename, "wb")
        tmpfile.write(response.content)
        tmpfile.close()
        shutil.copy(tmp_filename, "/tmp")
        tf = tarfile.open(tmp_filename)
        pl_path = os.path.join(output_directory, tf.next().name)
        tf.extractall(output_directory)
        os.unlink(tmp_filename)
        pangolearn_unpacked_dir = os.path.join(pl_path, "pangoLEARN")
        if v3datatree:
            # pangolin v3 expects a datadir with the entire pangoLEARN
            # subfolder in it.
            # In addition, it will only use the data if the __init__.py file
            # contained within that subfolder declares a __version__ that is
            # newer than the version installed with pangolin.
            pangolearn_dir = os.path.join(
                output_directory, tmp_filename, "pangoLEARN"
            )
            os.mkdir(os.path.dirname(pangolearn_dir))
            # rewrite the __init__.py file and make the __version__ string
            # appear newer than anything that might come with pangolin by
            # prepending a "v" to it.
            pangolearn_init = open(
                os.path.join(pangolearn_unpacked_dir, "__init__.py")
            ).readlines()
            with open(
                os.path.join(pangolearn_unpacked_dir, "__init__.py"), "w"
            ) as o:
                for line in pangolearn_init:
                    if line.startswith('__version__ = "'):
                        line = line.replace(
                            '__version__ = "', '__version__ = "v'
                        )
                    o.write(line)
        else:
            # Earlier versions of pangolin expect a datadir with just the
            # contents of the downloaded pangoLEARN subfolder in it and don't
            # care about the declared version in __init__.py.
            pangolearn_dir = os.path.join(
                output_directory, tmp_filename
            )
        os.rename(
            pangolearn_unpacked_dir,
            pangolearn_dir
        )
        shutil.rmtree(pl_path)
        return tmp_filename
    else:
        response.raise_for_status()


def parse_date(d):
    # Tries to parse the first 10 chars of d as a date, which currently
    # succeeds for all pangolearn model releases.
    return datetime.datetime.strptime(d[:10], "%Y-%m-%d")


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
        releases = filter_by_date([], start_date=args.start_date, end_date=args.end_date)
        for release in releases:
            print(release["tag_name"], release["tarball_url"].split("/")[-1], release["date"])
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
        releases = [next(get_model_list(existing_release_tags))]
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
        if args.pangolearn_format_version is not None:
            version = args.pangolearn_format_version
        else:
            # 2021-05-27 was the first release of pangoLEARN for pangolin 3, which changed DB format
            if release["date"] >= datetime.datetime(2021, 5, 27):
                version = '3.0'
            else:
                version = '1.0'
        fname = download_and_unpack(
            release["tarball_url"],
            output_directory,
            v3datatree=version == '3.0'
        )
        data_manager_dict["data_tables"][args.datatable_name].append(
            dict(
                value=release["tag_name"],
                description=release["name"],
                format_version=version,
                path=output_directory + "/" + fname,
            )
        )
    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
