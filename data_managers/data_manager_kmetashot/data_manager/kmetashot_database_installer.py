#!/usr/bin/env python

import argparse
import json
import os
import sys
from urllib.parse import urlparse
from urllib.request import HTTPError, Request, urlopen

#urls for from where to downlaod ref data
urls = {
    '1':{
        'path': "http://srv00.recas.ba.infn.it/webshare/brunofosso/kMetaShot_reference.h5",
        'release_date': "2022-07-31"
    },
    '2':{
        'path': "https://zenodo.org/records/17375120/files/kMetaShot_bacteria_archaea_2025-05-22.h5",
        'release_date': "2025-05-22"
    },
}


def is_urlfile(url):
    #check if file exist
    try:
        r = urlopen(url)  # response
        return r.getcode() < 400
    except HTTPError:
        return False


def url_download(url, target_directory):

    # download the url
    url_parts = urlparse(url)
    tarball = os.path.abspath(
        os.path.join(target_directory, os.path.basename(url_parts.path))
    )
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(tarball, "wb") as dst:
            while True:
                chunk = src.read(2**16)  # Read in 64 KB chunks instead of 1 KB
                if chunk:
                    dst.write(chunk)
                else:
                    break
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src is not None:
            src.close()

    return target_directory


def create_data_manager_entry(database_name, release, file_path):
    data_manager_entry = {}
    data_manager_entry["value"] = ( urls['release']['release_date']
    )
    data_manager_entry["name"] = f"kMetaShot reference data {urls['release']['release_date']}"
    data_manager_entry["path"] = file_path
    data_manager_entry["version"] = release
    return data_manager_entry


def download(release, test, out_file):

    with open(out_file) as fh:
        params = json.load(fh)

    target_directory = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_directory)

    if test:
        # make use of the test to check if all urls exists
        for _version, items in urls.items():
            for url in items.values():
                assert is_urlfile(url)

    data_manager_json = {"data_tables": {}}

    url = urls[release]["path"]
    file_path = url_download(url, target_directory)
    data_manager_json["data_tables"]["kmetashot"] = [
        create_data_manager_entry("Full Database", release, file_path)
    ]

    # store in dedicated metadata table
    with open(out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--release", dest="release", help="Release of the kMetaShot ref data"
    )
    parser.add_argument("--out_file", dest="out_file", help="JSON output file")
    parser.add_argument(
        "--test",
        dest="test",
        action="store_true",
        help="Run test",
    )
    args = parser.parse_args()

    download(
        args.release,
        args.test,
        args.out_file,
    )