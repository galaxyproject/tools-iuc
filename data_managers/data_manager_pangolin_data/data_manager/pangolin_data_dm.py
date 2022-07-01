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
from typing import Generator, List

import requests
from packaging import version


def parse_date(d: str) -> datetime.datetime:
    # Parses the publication date from the GitHub API or user input into a datetime object.
    date = None
    try:
        date = datetime.datetime.strptime(d, '%Y-%m-%dT%H:%M:%SZ')
    except ValueError:
        date = datetime.datetime.strptime(d, "%Y-%m-%d")
    return date


def get_model_list(
    existing_release_tags: List[str],
    package: str
) -> Generator[dict, None, None]:
    page_num = 0
    while True:
        url = f"https://api.github.com/repos/cov-lineages/{package}/releases"
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
                    date=parse_date(e["published_at"]),
                    tarball_url=e["tarball_url"],
                )
        else:
            response.raise_for_status()


def filter_by_date(existing_release_tags: List[str], package_name: str,
                   start_date: datetime.datetime = None, end_date: datetime.datetime = None) -> List[dict]:
    ret = []
    for release in get_model_list(existing_release_tags, package_name):
        if start_date and release["date"] < start_date:
            break
        if not end_date or release["date"] <= end_date:
            ret.append(release)

    return ret


def filter_by_version(existing_release_tags: List[str],
                      package_name: str, start_version: str, end_version: str) -> List[dict]:
    ret = []
    for release in get_model_list(existing_release_tags, package_name):
        if start_version is not None and version.parse(release["tag_name"]) < version.parse(start_version):
            # we can stop looking because releases are sorted by version, from newest to oldest
            break
        if end_version is None or version.parse(release["tag_name"]) <= version.parse(end_version):
            ret.append(release)
    return ret


def download_and_unpack(dependency: str, release: str, output_directory: str) -> pathlib.Path:
    url = f"git+https://github.com/cov-lineages/{dependency}.git@{release}"
    dependency_package_name = dependency.replace('-', '_')
    print(dependency, dependency_package_name, "release", release, "output_directory", output_directory)
    output_path = pathlib.Path(output_directory) / dependency_package_name / release
    print("output path is", output_path, "release is", release)
    with tempfile.TemporaryDirectory() as tmpdir:
        pip_command = [sys.executable, '-m', 'pip', 'install', '--upgrade',
                       '--target', tmpdir, url]
        # output is saved in tmpdir/dependency, final output needs to be
        # in output_directory/dependency/release
        subprocess.run(pip_command,
                       check=True)
        shutil.move(str(pathlib.Path(tmpdir) / dependency_package_name), str(output_path))
    return output_path


def comma_split(args: str) -> List[str]:
    return args.split(",")


def git_lfs_install():
    """
    'git-lfs install' must be run after installing git-lfs and before cloning a repo
    that uses Git LFS. Code taken from pangolin repo.
    """
    try:
        subprocess.run(['git-lfs', 'install'],
                       check=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode('utf-8')
        sys.stderr.write(f"Error: {e}:\n{stderr}\n")
        sys.exit(-1)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--testmode", default=False, action="store_true")
    parser.add_argument("--latest", default=False, action="store_true")
    parser.add_argument("--start_date", type=parse_date)
    parser.add_argument("--end_date", type=parse_date)
    parser.add_argument("--start_version", type=str)
    parser.add_argument("--end_version", type=str)
    parser.add_argument("--overwrite", default=False, action="store_true")
    parser.add_argument('--known_revisions', type=comma_split)
    parser.add_argument("datatable_name")
    parser.add_argument("datatable_cache_filename")
    parser.add_argument("galaxy_config")
    args = parser.parse_args()

    if args.datatable_name == 'pangolin_data':
        package_name = 'pangolin-data'
        min_version_key = 'min_pangolin_version'
        min_version = '4'
    elif args.datatable_name == 'pangolin_constellations':
        package_name = 'constellations'
        min_version_key = 'min_scorpio_version'
        min_version = '0'
    elif args.datatable_name == 'pangolin_assignment':
        package_name = 'pangolin-assignment'
        min_version_key = 'min_pangolin_version'
        min_version = '4'
        git_lfs_install()
    else:
        sys.exit(f"Unknown data table {args.datatable_name}")

    if args.testmode:
        if args.start_version is not None:
            releases = filter_by_version([], package_name, args.start_version, args.end_version)
        else:
            releases = filter_by_date([], package_name, start_date=args.start_date, end_date=args.end_date)
        for release in releases:
            print(release["tag_name"], release["tarball_url"].split("/")[-1], release["date"])
        sys.exit(0)

    with open(args.galaxy_config) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)

    try:
        with open(args.datatable_cache_filename) as fh:
            data_manager_dict = json.load(fh)
    except IOError:
        # on the first run this file doesn't exist
        data_manager_dict = {}

    if 'data_tables' in data_manager_dict:
        if args.datatable_name not in data_manager_dict['data_tables']:
            # got a data_tables entry, probably from a previous run of this script,
            # but no entry for this specific data table
            data_manager_dict['data_tables'][args.datatable_name] = []
    else:
        # got no entry for data tables, start from scratch
        data_manager_dict = {"data_tables": {args.datatable_name: []}}

    # known-revisions is populated from the Galaxy `pangolin_data` data table by the wrapper
    if args.known_revisions is not None:
        existing_release_tags = set(args.known_revisions)
    else:
        existing_release_tags = set()
    if args.latest:
        latest_release = next(get_model_list([], package_name))
        if latest_release["tag_name"] in existing_release_tags:
            releases = []
        else:
            releases = [latest_release]
    elif args.start_version is not None or args.end_version is not None:
        releases = filter_by_version(existing_release_tags,
                                     package_name,
                                     args.start_version, args.end_version)
    else:
        releases = filter_by_date(
            existing_release_tags,
            package_name,
            start_date=args.start_date, end_date=args.end_date
        )
    releases_to_download = [
        release
        for release in releases
        if release["tag_name"] not in existing_release_tags
    ]
    for release in releases_to_download:

        fname = download_and_unpack(
            package_name,
            release['tag_name'],
            output_directory
        )
        if fname is not None:
            data_manager_dict["data_tables"][args.datatable_name].append(
                {
                    'value': release["tag_name"],
                    'description': release["name"],
                    min_version_key: min_version,
                    'date': release["date"].isoformat(),  # ISO 8601 is easily sortable
                    'path': str(output_directory / fname)
                }
            )
    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.datatable_cache_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
