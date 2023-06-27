#!/usr/bin/env python

import argparse
import datetime
import json
import operator
import pathlib
import shutil
import subprocess
import sys
import tempfile
from io import StringIO
from typing import Generator, TextIO

import requests


def parse_date(d: str) -> datetime.datetime:
    # Parses the publication date from the GitHub API or user input into a datetime object.
    date = None
    try:
        date = datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        date = datetime.datetime.strptime(d, "%Y-%m-%d")
    return date


def get_model_list(package: str) -> Generator[dict, None, None]:
    page_num = 0
    while True:
        url = f"https://api.github.com/repos/cov-lineages/{package}/releases"
        page_num += 1
        response = requests.get(url + f"?page={page_num}")
        if response.status_code == 200:
            release_list_chunk = json.loads(response.text)
            if not release_list_chunk:
                # past the last page of results
                return
            for e in release_list_chunk:
                if e["prerelease"]:
                    continue
                yield dict(
                    tag_name=e["tag_name"],
                    name=e["name"],
                    date=parse_date(e["published_at"]),
                    tarball_url=e["tarball_url"],
                )
        else:
            response.raise_for_status()


def download_and_unpack(
    dependency: str, release: str, output_directory: str
) -> pathlib.Path:
    url = f"git+https://github.com/cov-lineages/{dependency}.git@{release}"
    dependency_package_name = dependency.replace("-", "_")
    output_path = pathlib.Path(output_directory) / dependency_package_name / release
    with tempfile.TemporaryDirectory() as tmpdir:
        pip_command = [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--upgrade",
            "--target",
            tmpdir,
            url,
        ]
        # output is saved in tmpdir/dependency, final output needs to be
        # in output_directory/dependency/release
        subprocess.run(pip_command, check=True)
        shutil.move(
            str(pathlib.Path(tmpdir) / dependency_package_name), str(output_path)
        )
    return output_path


def fetch_compatibility_info(
    package_name: str,
    url: str = "https://raw.githubusercontent.com/cov-lineages/pangolin/master/pangolin/data/data_compatibility.csv",
) -> list[dict[str, str]]:
    response = requests.get(url)
    if response.status_code == 200:
        compatibility = read_compatibility_info(StringIO(response.text), package_name)
        return compatibility
    else:
        return {}


def read_compatibility_info(
    input_file: TextIO, package_name: str
) -> list[dict[str, str]]:
    compatibility = {}
    for line in input_file:
        fields = line.strip().split(",")
        if fields[0] != package_name:
            continue
        if package_name == "constellations":
            compatibility[fields[1]] = fields[3]
        else:
            # for pangolin-data and pangolin-assignment
            compatibility[fields[1]] = fields[2]
    return compatibility


def comma_split(args: str) -> list[str]:
    return args.split(",")


def git_lfs_install():
    """
    'git-lfs install' must be run after installing git-lfs and before cloning a repo
    that uses Git LFS. Code taken from pangolin repo.
    """
    try:
        subprocess.run(
            ["git-lfs", "install"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode("utf-8")
        sys.stderr.write(f"Error: {e}:\n{stderr}\n")
        sys.exit(-1)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--latest", default=False, action="store_true")
    parser.add_argument("--version_compatibility_file", type=argparse.FileType())
    parser.add_argument("--versions", type=comma_split)
    parser.add_argument("--overwrite", default=False, action="store_true")
    parser.add_argument("--known_revisions", type=comma_split)
    parser.add_argument("datatable_name")
    parser.add_argument("datatable_cache_filename")
    parser.add_argument("galaxy_config")
    args = parser.parse_args()

    if args.datatable_name == "pangolin_data":
        package_name = "pangolin-data"
        min_version_key = "min_pangolin_version"
    elif args.datatable_name == "pangolin_constellations":
        package_name = "constellations"
        min_version_key = "min_scorpio_version"
    elif args.datatable_name == "pangolin_assignment":
        package_name = "pangolin-assignment"
        min_version_key = "min_pangolin_version"
        git_lfs_install()
    else:
        sys.exit(f"Unknown data table {args.datatable_name}")

    with open(args.galaxy_config) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)

    try:
        with open(args.datatable_cache_filename) as fh:
            data_manager_dict = json.load(fh)
    except IOError:
        # on the first run this file doesn't exist
        data_manager_dict = {}

    if "data_tables" in data_manager_dict:
        if args.datatable_name not in data_manager_dict["data_tables"]:
            # got a data_tables entry, probably from a previous run of this script,
            # but no entry for this specific data table
            data_manager_dict["data_tables"][args.datatable_name] = []
    else:
        # got no entry for data tables, start from scratch
        data_manager_dict = {"data_tables": {args.datatable_name: []}}

    # known-revisions is populated from the Galaxy `pangolin_data` data table by the wrapper
    if args.known_revisions is not None:
        existing_release_tags = set(args.known_revisions)
    else:
        existing_release_tags = set()
    if args.latest:
        compatibility = fetch_compatibility_info(package_name)
        for latest_release in get_model_list(package_name):
            # choose the first release for which we have compatibility info
            version = latest_release["tag_name"].lstrip("v.")
            if version in compatibility:
                latest_release[min_version_key] = compatibility[version]
                break
        if latest_release["tag_name"] in existing_release_tags:
            releases = []
        else:
            releases = [latest_release]
    else:
        compatibility = read_compatibility_info(
            args.version_compatibility_file, package_name
        )
        downloadable_releases = get_model_list(package_name)
        releases_wanted = set(args.versions) - set(
            [tag.lstrip("v.") for tag in existing_release_tags]
        )
        releases = []
        for release in downloadable_releases:
            version = release["tag_name"].lstrip("v.")
            if version in releases_wanted:
                if version in compatibility:
                    # only add the releases for which we have compatibility info
                    release[min_version_key] = compatibility[version]
                    releases.append(release)
                    releases_wanted.remove(version)
                    if not releases_wanted:
                        # we've found all the releases we want
                        break
        if releases_wanted:
            missing_releases = " ".join(releases_wanted)
            sys.exit(
                f"Some of the requested releases ({missing_releases}) are not available."
            )

    for release in releases:
        fname = download_and_unpack(package_name, release["tag_name"], output_directory)
        if fname is not None:
            data_manager_dict["data_tables"][args.datatable_name].append(
                {
                    "value": release["tag_name"],
                    "description": release["name"],
                    min_version_key: release[min_version_key],
                    "date": release["date"].isoformat(),  # ISO 8601 is easily sortable
                    "path": str(output_directory / fname),
                }
            )
    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.datatable_cache_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
