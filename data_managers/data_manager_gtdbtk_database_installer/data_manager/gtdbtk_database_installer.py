#!/usr/bin/env python

import argparse
import gzip
import json
import os
import shutil
import sys
import tarfile
from datetime import datetime
from urllib.parse import urlparse
from urllib.request import Request, urlopen

# rather provide the urls based on the release, less error potential for the admins !
urls = {
    "202": {
        "full": "https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz",
        "meta_ar": "https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz",
        "meta_bac": "https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz",
    },
    "207": {
        "full": "https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_data.tar.gz",
        "meta_ar": "https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz",
        "meta_bac": "https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz",
    },
    "214": {
        "full": "https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz",
        "meta_ar": "https://data.gtdb.ecogenomic.org/releases/release214/214.1/ar53_metadata_r214.tsv.gz",
        "meta_bac": "https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_metadata_r214.tsv.gz",
    },
    "220": {
        "full": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz",
        "meta_ar": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/ar53_metadata_r220.tsv.gz",
        "meta_bac": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz",
    },
    "test": {  # using VERSION to check if files are there
        "full": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/VERSION.txt",
        "meta_ar": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/ar53_metadata_r220.tsv.gz",
        "meta_bac": "https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz",
    },
}


def url_download(url, target_directory):
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
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src is not None:
            src.close()
    if tarfile.is_tarfile(tarball):
        fh = tarfile.open(tarball, "r:*")
    else:
        # unzip metadata file
        if ".gz" in tarball:
            with gzip.open(tarball, "rb") as f_in:
                unzipped_file = tarball.strip(".gz")
                with open(unzipped_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(tarball)
                folder_of_unzipped_file = os.path.dirname(unzipped_file)
            return folder_of_unzipped_file
        else:
            # this is basically only the return for the test not using a tarfile
            return tarball
    fh.extractall(target_directory)
    fh.close()
    os.remove(tarball)
    # The tarball extraction will create a directory named
    # something like release202 in the target_directory, so
    # we need to move the items in that directory to the
    # target directory.
    subdir = next(os.walk(target_directory))[1][0]
    subdir_path = os.path.join(target_directory, subdir)
    items = os.listdir(subdir_path)
    for item in items:
        item_path = os.path.join(subdir_path, item)
        shutil.move(item_path, target_directory)
    os.rmdir(subdir_path)
    return target_directory


def download(database_name, release, meta, test, out_file):

    with open(out_file) as fh:
        params = json.load(fh)

    target_directory = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_directory)

    if test:
        release = "test"

    # download both taxonomy metadata tables
    if meta:
        url = urls[release]["meta_ar"]
        file_path = url_download(url, target_directory)
        url = urls[release]["meta_bac"]
        file_path = url_download(url, target_directory)
    # download the full DB
    else:
        url = urls[release]["full"]
        file_path = url_download(url, target_directory)

    time = datetime.utcnow().strftime("%Y-%m-%d")

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry["value"] = f"{database_name}_release_{release}_downloaded_{time}"
    data_manager_entry["name"] = database_name
    data_manager_entry["path"] = file_path
    data_manager_entry["version"] = release

    # store in dedicated metadata table
    if meta:
        data_manager_json["data_tables"][
            "gtdbtk_database_metadata_versioned"
        ] = data_manager_entry
    else:
        data_manager_json["data_tables"][
            "gtdbtk_database_versioned"
        ] = data_manager_entry

    with open(out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--database_name", dest="database_name", help="GTDB-Tk database display name"
)

parser.add_argument("--version", dest="version", help="DB version")

parser.add_argument(
    "--release", dest="release", help="Release of the GTDB-Tk database version"
)
parser.add_argument("--out_file", dest="out_file", help="JSON output file")
parser.add_argument(
    "--meta",
    dest="meta",
    action="store_true",
    help="Store meta data flag",
)

parser.add_argument(
    "--test",
    dest="test",
    action="store_true",
    help="Run test",
)

args = parser.parse_args()

download(
    args.database_name,
    args.release,
    args.meta,
    args.test,
    args.out_file,
)
