#!/usr/bin/env python

import argparse
import gzip
import hashlib
import json
import os
import posixpath
import shutil
import sys
import tarfile
from urllib.parse import urlparse
from urllib.request import HTTPError, Request, urlopen

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
    "226": {
        "full": "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz",
        "meta_ar": "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_metadata_r226.tsv.gz",
        "meta_bac": "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_metadata_r226.tsv.gz",
    },
    "232.0": {
        "full": "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz",
        "meta_ar": "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/ar53_metadata_r232.tsv.gz",
        "meta_bac": "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/bac120_metadata_r232.tsv.gz",
    },
    "mocked_232": {
        "full": "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/mockup_db/r232_mockup.tar.gz",
    },
}

# Name of the checksum files published by GTDB in the release directories
CHECKSUM_FILE_NAMES = ("MD5SUM.txt", "MD5SUM")
# Number of directory levels to search upwards for a checksum file
CHECKSUM_SEARCH_LEVELS = 5

# Test fixtures: a miniature database (1.7 kB), one without and one with a
# published MD5 sum (https://zenodo.org/records/23918153/files/MD5SUM.txt).
# In test mode release 232 uses the fixture with MD5 sum, all other releases the
# one without, such that both cases are covered by the tests.
TEST_DB_URL = "https://zenodo.org/records/13734217/files/release220-test.tar.gz"
TEST_DB_URL_MD5 = "https://zenodo.org/records/23918153/files/release220-test.tar.gz"
TEST_DB_URLS = {"232.0": TEST_DB_URL_MD5}


def is_urlfile(url):
    # Check if online file exists
    try:
        r = urlopen(url)  # response
        try:
            return r.getcode() < 400
        finally:
            r.close()
    except HTTPError:
        return False


def find_checksum_url(url, cache):
    """Return the URL of the MD5 checksum file that covers ``url`` (or None).

    GTDB publishes MD5 sums in the root directory of each release, e.g.
    .../releases/release232/232.0/MD5SUM.txt, while the files themselves may
    live in sub directories. Hence, we probe ``url``'s directory and its
    parents. Results are cached per directory to avoid repeated requests.
    """
    url_parts = urlparse(url)
    directory = posixpath.dirname(url_parts.path)
    for _ in range(CHECKSUM_SEARCH_LEVELS):
        if directory in cache:
            found = cache[directory]
        else:
            found = None
            for filename in CHECKSUM_FILE_NAMES:
                candidate = url_parts._replace(
                    path=posixpath.join(directory, filename)
                ).geturl()
                if is_urlfile(candidate):
                    found = candidate
                    break
            cache[directory] = found
        if found:
            return found
        parent = posixpath.dirname(directory)
        if parent == directory:
            break
        directory = parent
    return None


def fetch_md5sums(checksum_url):
    """Return a {relative file path: md5} mapping from a checksum file."""
    md5sums = {}
    with urlopen(Request(checksum_url)) as fh:
        content = fh.read().decode("utf-8", errors="replace")
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        md5sum, filename = fields[0], fields[-1]
        if filename.startswith("*"):  # binary mode marker written by md5sum
            filename = filename[1:]
        if filename.startswith("./"):
            filename = filename[2:]
        md5sums[filename] = md5sum.lower()
    return md5sums


def expected_md5(url, checksum_cache, md5sums_cache):
    """Return ``(relative path, md5, checksum url)`` published for ``url``.

    Files for which the data source does not publish a checksum (e.g. the
    mockup database used for testing) yield None.
    """
    checksum_url = find_checksum_url(url, checksum_cache)
    if checksum_url is None:
        return None
    if checksum_url not in md5sums_cache:
        try:
            md5sums_cache[checksum_url] = fetch_md5sums(checksum_url)
        except Exception as e:  # checksums are a bonus, do not fail on them
            print(f"Could not read {checksum_url}: {e}", file=sys.stderr)
            md5sums_cache[checksum_url] = {}
    md5sums = md5sums_cache[checksum_url]

    path = urlparse(url).path
    base = posixpath.dirname(urlparse(checksum_url).path)
    # checksum files list paths relative to their own directory
    if path.startswith(base + "/"):
        relpath = path[len(base) + 1:]
        if relpath in md5sums:
            return relpath, md5sums[relpath], checksum_url
    # fall back to a unique match of the file name
    basename = posixpath.basename(path)
    candidates = {
        relpath: md5sum
        for relpath, md5sum in md5sums.items()
        if posixpath.basename(relpath) == basename
    }
    if len(set(candidates.values())) == 1:
        relpath, md5sum = next(iter(candidates.items()))
        return relpath, md5sum, checksum_url
    return None


def extract_tar_iteratively(tarball, target_directory):
    """
    Extracts a .tar, .tar.gz, or .tar.bz2 archive iteratively in a memory-efficient manner.

    This function processes the contents of the archive member-by-member, ensuring only
    one file or directory is loaded into memory at any given time. It handles the creation
    of directories and symbolic links, and streams large files to disk in chunks to avoid
    memory overload.

    Args:
        tarball (str): Path to the tar archive (e.g., .tar, .tar.gz, .tar.bz2) to be extracted.
        target_directory (str): The destination directory where the archive content
                                will be extracted.

    Raises:
        OSError: If there is an issue with file or directory creation, or writing to disk.
        tarfile.TarError: If there is an issue opening or reading the tar archive.

    Example Usage:
        extract_tar_iteratively("archive.tar.gz", "/path/to/extract")

    Notes:
        - The function supports symbolic and hard links present in the tar archive.
        - It ensures that directories are created before files are extracted.
        - Large files are streamed to disk in 1 MB chunks to minimize memory usage.
        - This function does not return anything but will populate the target directory with
          the extracted content.
    """

    with tarfile.open(tarball, "r:*") as fh:
        for member in fh:
            # Full path to where the member should be extracted
            member_path = os.path.join(target_directory, member.name)

            if member.isdir():
                # If it's a directory, ensure it exists
                os.makedirs(member_path, exist_ok=True)
            elif member.isfile():
                # If it's a file, extract it in chunks to avoid memory spikes
                with fh.extractfile(member) as source, open(
                    member_path, "wb"
                ) as target:
                    shutil.copyfileobj(
                        source, target, length=1024 * 1024
                    )  # 1 MB chunks
            elif member.issym() or member.islnk():
                # Handle symlinks or hard links if necessary
                target_link = os.path.join(target_directory, member.name)
                if member.issym():
                    os.symlink(member.linkname, target_link)
                elif member.islnk():
                    os.link(member.linkname, target_link)


def url_download(url, target_directory, meta, md5sum=None):

    # download the url
    url_parts = urlparse(url)
    tarball = os.path.abspath(
        os.path.join(target_directory, os.path.basename(url_parts.path))
    )
    src = None
    dst = None
    md5 = hashlib.md5()
    try:
        req = Request(url)
        src = urlopen(req)
        with open(tarball, "wb") as dst:
            while True:
                chunk = src.read(2**16)  # Read in 64 KB chunks instead of 1 KB
                if chunk:
                    md5.update(chunk)
                    dst.write(chunk)
                else:
                    break
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src is not None:
            src.close()

    # verify the download against the checksum published by GTDB
    if md5sum is not None and md5.hexdigest() != md5sum:
        os.remove(tarball)
        sys.exit(
            f"Checksum mismatch for {url}: expected {md5sum} but found "
            f"{md5.hexdigest()}. The file has not been installed."
        )

    # extract the metadata
    if meta:
        # extract the content of *.tar.gz into the target dir
        if tarfile.is_tarfile(tarball):
            extract_tar_iteratively(tarball, target_directory)
            os.remove(tarball)
            return target_directory  # return path to output folder
        # extract the content of *.gz into the target dir
        elif ".gz" in tarball:
            with gzip.open(tarball, "rb") as f_in:
                unzipped_file = tarball.strip(".gz")
                with open(unzipped_file, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    os.remove(tarball)
                    folder_of_unzipped_file = os.path.dirname(unzipped_file)
            return folder_of_unzipped_file
        else:
            sys.exit(
                "No correct input format for metadata file, must be .tar.gz or .gz"
            )
    else:
        # handle the DB
        # extract the content of the folder in the tar.gz into the target dir
        if tarfile.is_tarfile(tarball):
            extract_tar_iteratively(tarball, target_directory)
            os.remove(tarball)
        else:
            # handle the test case for the DB
            return tarball

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


def version_of(release):
    """Return the version of a release, i.e. the release without minor version.

    The keys of the ``urls`` table are the releases as GTDB names their download
    directories (e.g. ``232.0``). This is used for the value and the name of the
    data table entry, but the GTDB-Tk tools filter on the plain release (e.g.
    ``232``), which GTDB keeps stable when it republishes the files of a release.
    """
    return release.split(".")[0]


def create_data_manager_entry(database_name, release, file_path):
    data_manager_entry = {}
    data_manager_entry["value"] = (
        f"{database_name.replace(' ', '_').lower()}_release_{release}"
    )
    if release == "mocked_232":
        data_manager_entry["name"] = "Mocked GTBD DB (232)"
    else:
        data_manager_entry["name"] = f"{database_name} - release {release}"
    data_manager_entry["path"] = file_path
    data_manager_entry["version"] = version_of(release)
    return data_manager_entry


def download(release, meta, test, out_file):

    with open(out_file) as fh:
        params = json.load(fh)

    target_directory = params["output_data"][0]["extra_files_path"]
    os.makedirs(target_directory)

    if test:
        # switch the DB to use one of the test cases
        urls[release]["full"] = TEST_DB_URLS.get(release, TEST_DB_URL)

        # make use of the test to check if all urls exists
        for _version, items in urls.items():
            for url in items.values():
                assert is_urlfile(url)

    data_manager_json = {"data_tables": {}}

    checksum_cache = {}
    md5sums_cache = {}

    def download_verified(url):
        """Download url, verifying it against the checksum of the data source.

        Reports which checksum was used, and returns the local path of the
        installed file.
        """
        filename = posixpath.basename(urlparse(url).path)
        expected = expected_md5(url, checksum_cache, md5sums_cache)
        if expected is None:
            print(
                f"No remote MD5 sum available for {filename}: "
                "download not verified"
            )
            return url_download(url, target_directory, meta)
        _, md5sum, checksum_url = expected
        file_path = url_download(url, target_directory, meta, md5sum=md5sum)
        print(f"Remote MD5 sum used for {filename}: {md5sum} (from {checksum_url})")
        return file_path

    # download taxonomy metadata tables
    if meta:
        download_verified(urls[release]["meta_ar"])
        file_path = download_verified(urls[release]["meta_bac"])

        data_manager_json["data_tables"]["gtdbtk_database_metadata_versioned"] = [
            create_data_manager_entry("Metadata Tables", release, file_path)
        ]
    # download the full DB
    else:
        file_path = download_verified(urls[release]["full"])
        data_manager_json["data_tables"]["gtdbtk_database_versioned"] = [
            create_data_manager_entry("Full Database", release, file_path)
        ]

    # store in dedicated metadata table
    with open(out_file, "w") as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

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
        args.release,
        args.meta,
        args.test,
        args.out_file,
    )
