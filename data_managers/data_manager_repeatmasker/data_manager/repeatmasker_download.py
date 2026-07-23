#!/usr/bin/env python
"""Download a Dfam FamDB library for RepeatMasker and register it in a data table.

Dfam distributes its libraries in the partitioned FamDB (v3) format: a single
required *root* partition plus one file per *component* partition.  The
components are split by curation status (curated/uncurated) and model type
(consensus sequences, used by the rmblast search engine, or profile HMMs, used
by nhmmer).  This data manager downloads the root partition (always) together
with whichever components the admin selects, laying the ``*.h5`` files out in a
single directory that RepeatMasker (>= 4.2.4) can consume via ``-libdir``.

See https://www.dfam.org/releases/current/families/FamDB/README.txt
"""

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
from datetime import date
from urllib.error import HTTPError
from urllib.request import Request, urlopen

# Map a Dfam release to its release directory and FamDB file-name prefix.
# Adding a new release is a one-line change here (the file layout is stable).
RELEASES = {
    "4.0": {"dir": "Dfam_4.0", "prefix": "dfam40"},
    "3.9": {"dir": "Dfam_3.9", "prefix": "dfam39"},
    "3.8": {"dir": "Dfam_3.8", "prefix": "dfam38"},
}

BASE_URL = "https://www.dfam.org/releases/{dir}/families/FamDB/"

# component key -> FamDB file-name infix
COMPONENTS = {
    "curated_consensus": "curated.consensus",
    "uncurated_consensus": "uncurated.consensus",
    "curated_hmm": "curated.hmm",
    "uncurated_hmm": "uncurated.hmm",
}

# A tiny (~0.3 MB) real partition used to exercise the full download / gunzip /
# checksum / registration path quickly during automated testing.
TEST_COMPONENT = "uncurated_consensus"
CHUNK = 2 ** 16


def url_exists(url):
    """Return True if a HEAD-like GET on url succeeds with status < 400."""
    try:
        return urlopen(Request(url)).getcode() < 400
    except HTTPError:
        return False


def fetch_md5(url):
    """Return the expected md5 hex digest from Dfam's ``<file>.md5`` sidecar."""
    with urlopen(Request(url + ".md5")) as response:
        # sidecar format: "<md5>  <filename>"
        return response.read().decode().split()[0]


def download_and_extract(url, target_directory):
    """Download a gzipped ``.h5.gz`` partition, verify its md5, gunzip it.

    The uncompressed ``.h5`` file is written into ``target_directory`` and the
    compressed download is removed.  Raises on a checksum mismatch so a
    truncated or corrupted download never registers as a usable library.
    """
    gz_name = url.rsplit("/", 1)[1]
    gz_path = os.path.join(target_directory, gz_name)
    h5_path = os.path.join(target_directory, gz_name[: -len(".gz")])

    expected_md5 = fetch_md5(url)
    md5 = hashlib.md5()
    with urlopen(Request(url)) as src, open(gz_path, "wb") as dst:
        while True:
            chunk = src.read(CHUNK)
            if not chunk:
                break
            md5.update(chunk)
            dst.write(chunk)
    if md5.hexdigest() != expected_md5:
        sys.exit(
            "Checksum mismatch for {}: expected {}, got {}".format(
                gz_name, expected_md5, md5.hexdigest()
            )
        )

    with gzip.open(gz_path, "rb") as f_in, open(h5_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out, CHUNK)
    os.remove(gz_path)
    return os.path.basename(h5_path)


def iter_partition_urls(base_url, prefix, infix):
    """Yield partition URLs for a component, probing 0,1,2,... until one is absent.

    Dfam numbers component partitions contiguously from 0, so the first missing
    index marks the end.  This adapts automatically to how many partitions a
    component has in a given release instead of hard-coding per-release counts.
    """
    part = 0
    while True:
        url = "{}{}.{}.{}.h5.gz".format(base_url, prefix, infix, part)
        if not url_exists(url):
            break
        yield url
        part += 1


def download(release, components, test, out_file):
    if release not in RELEASES:
        sys.exit(
            "Unknown Dfam release '{}'. Known releases: {}".format(
                release, ", ".join(sorted(RELEASES))
            )
        )
    rel = RELEASES[release]
    base_url = BASE_URL.format(dir=rel["dir"])
    prefix = rel["prefix"]

    with open(out_file) as fh:
        params = json.load(fh)
    target_directory = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_directory)

    if test:
        # Exercise the real code path cheaply: fetch only the smallest partition
        # of a single component, and skip the large root/curated downloads.
        components = [TEST_COMPONENT]
        infix = COMPONENTS[TEST_COMPONENT]
        installed = [
            download_and_extract(
                "{}{}.{}.0.h5.gz".format(base_url, prefix, infix), target_directory
            )
        ]
    else:
        # The root partition is always required.
        installed = [
            download_and_extract(
                "{}{}.0.h5.gz".format(base_url, prefix), target_directory
            )
        ]
        for component in components:
            infix = COMPONENTS[component]
            for url in iter_partition_urls(base_url, prefix, infix):
                installed.append(download_and_extract(url, target_directory))

    today = date.today().strftime("%Y-%m-%d")
    selected = "+".join(components) if components else "root"
    entry = {
        "value": "dfam_{}_{}_{}".format(release, selected, today).replace(".", "_"),
        "name": "Dfam {} ({}) [{}]".format(release, ", ".join(components) or "root only", today),
        "version": release,
        "path": target_directory,
    }
    data_manager_json = {"data_tables": {"repeatmasker_famdb": [entry]}}

    sys.stderr.write(
        "Installed {} FamDB partition file(s) into {}\n".format(
            len(installed), target_directory
        )
    )
    with open(out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", required=True, help="Dfam release, e.g. 4.0")
    parser.add_argument(
        "--components",
        default="",
        help="comma-separated component keys to install in addition to the root partition",
    )
    parser.add_argument("--out_file", required=True, help="JSON output file")
    parser.add_argument(
        "--test", action="store_true", help="download only a tiny partition for testing"
    )
    args = parser.parse_args()

    components = [c for c in args.components.split(",") if c]
    unknown = [c for c in components if c not in COMPONENTS]
    if unknown:
        sys.exit("Unknown component(s): {}".format(", ".join(unknown)))

    download(args.release, components, args.test, args.out_file)
