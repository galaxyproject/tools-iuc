#!/usr/bin/env python

import argparse
import hashlib
import json
import operator
import os
import re
import shutil
import subprocess
import sys
import tarfile

import requests


GH_REPO_API = 'https://api.github.com/repos/ebi-pf-team/interproscan/'
MD5_URL = 'http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{version}/interproscan-{version}-64-bit.tar.gz.md5'
DATA_URL = 'http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{version}/interproscan-{version}-64-bit.tar.gz'

# For tests: download a smaller archive containing *some* data
PARTIAL_URL = 'https://github.com/ebi-pf-team/interproscan/archive/{version}.tar.gz'


def list_tags(url=None):

    if not url:
        url = GH_REPO_API + 'tags'

    resp = requests.get(url=url)
    data = resp.json()

    tags = []
    for tag in data:
        if re.match(r"^[0-9]\.[0-9]{2}-[0-9]{2}\.[0-9]$", tag['name']):
            tags.append(tag['name'])

    if 'next' in resp.links:
        tags += list_tags(resp.links['next']['url'])

    return sorted(tags)


def download_file(url, dest):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def main():
    parser = argparse.ArgumentParser(description='Download data for InterProScan')
    parser.add_argument('--partial', dest='partial', action='store_true', help='Only download a small subset of data (for testing)')
    parser.add_argument('-v', '--version', help='Specify an InterProScan version (default: latest)')
    parser.add_argument("datatable_name")
    parser.add_argument("galaxy_datamanager_filename")

    args = parser.parse_args()

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][args.datatable_name] = data_manager_dict[
        "data_tables"
    ].get(args.datatable_name, [])

    os.mkdir(output_directory)

    all_tags = list_tags()

    if args.version:
        if args.version not in all_tags:
            raise RuntimeError("Version '%s' is not valid" % args.version)
        tag = args.version
    else:
        tag = all_tags[-1]

    setup_script = 'initial_setup.py'
    sub_version = re.match(r"^[0-9]\.([0-9]{2})-[0-9]{2}\.[0-9]$", tag)
    if sub_version and len(sub_version.groups()) == 1 and int(sub_version.group(1)) >= 58:
        # The setup script was renamed in 5.58
        setup_script = 'setup.py'
    else:
        raise RuntimeError("Sorry, this data manager can only download data for InterProScan >= 5.58-91.0. Use the 0.0.2 version for older versions of InterProScan.")

    print("Will download data for InterProScan version: %s" % tag)

    print("Getting MD5 checksum:")
    md5 = requests.get(url=MD5_URL.format(version=tag)).text
    if not re.match(r"^([a-fA-F\d]{32})  interproscan-[0-9]\.[0-9]{2}-[0-9]{2}\.[0-9]-64-bit.tar.gz$", md5):
        raise RuntimeError("Got invalid MD5 from the InterProScan FTP server: '%s'" % md5)
    print("%s" % md5)

    if args.partial:
        print("Downloading partial data tarball...")
        dest_tar = os.path.join(output_directory, PARTIAL_URL.format(version=tag).split('/')[-1])
        download_file(PARTIAL_URL.format(version=tag), dest_tar)
    else:
        print("Downloading data tarball...")
        dest_tar = os.path.join(output_directory, DATA_URL.format(version=tag).split('/')[-1])
        download_file(DATA_URL.format(version=tag), dest_tar)

        print("Finished, now checking md5...")
        m = hashlib.md5()
        blocksize = 2**20
        with open(dest_tar, 'rb') as tarball:
            while True:
                buf = tarball.read(blocksize)
                if not buf:
                    break
                m.update(buf)
        md5_computed = m.hexdigest()
        if not md5.startswith(md5_computed):
            raise RuntimeError("MD5 check failed: computed '%s', expected '%s'" % (md5_computed, md5))

    print("Ok, now extracting data...")
    tar = tarfile.open(dest_tar, "r:gz")
    tar.extractall(output_directory)
    tar.close()

    if args.partial:
        print("Moving partial data files around...")
        shutil.move(os.path.join(output_directory, 'interproscan-%s' % tag, 'core/jms-implementation/support-mini-x86-32/data/'), os.path.join(output_directory, 'data'))
    else:
        print("Moving data files around...")
        shutil.move(os.path.join(output_directory, 'interproscan-%s' % tag, 'data'), os.path.join(output_directory, 'data'))

    print("Done, removing tarball and unneeded files...")
    os.remove(dest_tar)
    shutil.rmtree(os.path.join(output_directory, 'interproscan-%s' % tag))

    print("Running {} (index hmm models)...".format(setup_script))
    # Write a temp properties file in work dir
    prop_file_src = os.path.join(os.path.dirname(os.path.realpath(shutil.which("interproscan.sh"))), 'interproscan.properties')
    with open(prop_file_src, 'r') as prop:
        prop_content = prop.read()
    prop_content = re.sub(r'^data\.directory=.*$', 'data.directory=%s' % os.path.join(output_directory, 'data'), prop_content, flags=re.M)
    with open('interproscan.properties', 'w') as prop:
        prop.write(prop_content)
    # Run the index command
    cmd_args = [os.path.join(os.path.dirname(os.path.realpath(shutil.which("interproscan.sh"))), setup_script), 'interproscan.properties']
    proc = subprocess.Popen(args=cmd_args, shell=False)
    out, err = proc.communicate()
    print(out)
    print(err, file=sys.stderr)
    return_code = proc.wait()
    if return_code:
        print("Error running {}.".format(setup_script), file=sys.stderr)
        sys.exit(return_code)

    data_manager_dict["data_tables"][args.datatable_name].append(
        dict(
            value=tag,
            description="InterProScan %s" % tag,
            interproscan_version=tag,
            path=output_directory,
        )
    )

    print("Saving data table content...")

    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)

    print("Finished.")


if __name__ == "__main__":
    main()
