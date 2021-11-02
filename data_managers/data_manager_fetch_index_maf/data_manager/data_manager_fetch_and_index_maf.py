#!/usr/bin/env python
import bz2
import ftplib
import gzip
import json
import optparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
import zipfile
from binascii import hexlify

CHUNK_SIZE = 2**20

DEFAULT_DATA_TABLE_NAME = "indexed_maf_files"

# Nice solution to opening compressed files (zip/bz2/gz) transparently
# https://stackoverflow.com/a/13045892/638445


class CompressedFile(object):
    magic = None
    file_type = None
    mime_type = None
    proper_extension = None

    def __init__(self, f):
        # f is an open file or file like object
        self.f = f
        self.accessor = self.open()

    @classmethod
    def is_magic(self, data):
        return hexlify(data).startswith(hexlify(self.magic))

    def open(self):
        return None


class ZIPFile(CompressedFile):
    magic = b'\x50\x4b\x03\x04'
    file_type = 'zip'
    mime_type = 'compressed/zip'

    def open(self):
        return zipfile.ZipFile(self.f)


class BZ2File(CompressedFile):
    magic = b'\x42\x5a\x68'
    file_type = 'bz2'
    mime_type = 'compressed/bz2'

    def open(self):
        return bz2.BZ2File(self.f)


class GZFile(CompressedFile):
    magic = b'\x1f\x8b\x08'
    file_type = 'gz'
    mime_type = 'compressed/gz'

    def open(self):
        return gzip.GzipFile(self.f)


# Factory function to create a suitable instance for accessing files
def get_compressed_file(filename):
    with open(filename, 'rb') as f:
        start_of_file = f.read(16)
        f.seek(0)
        for cls in (ZIPFile, BZ2File, GZFile):
            if cls.is_magic(start_of_file):
                f.close()
                return cls(filename)

        return None


def url_download(url, tmp=False, localpath=None):
    """Attempt to download file from a given url
    :param url: full url to file
    :type url: str.
    :returns: name of downloaded file
    :raises: ContentDecodingError, IOError
    """

    # Generate file_name
    file_name = url.split('/')[-1]
    if tmp:
        file_name = os.path.join(tempfile.mkdtemp(), file_name)
    elif localpath is not None:
        file_name = os.path.join(localpath, file_name)

    try:
        # download URL (FTP and HTTP work, probably local and data too)
        urllib.request.urlretrieve(url, file_name)

        # uncompress file if needed
        cf = get_compressed_file(file_name)
        if cf is not None:
            uncompressed_file_name = os.path.splitext(file_name)[0]
            with open(uncompressed_file_name, 'wb') as uncompressed_file:
                shutil.copyfileobj(cf.accessor, uncompressed_file)
            os.remove(file_name)
            file_name = uncompressed_file_name
    except IOError as e:
        sys.stderr.write('Error occured downloading reference file: %s' % e)
        os.remove(file_name)
    return file_name


def generate_metadata(params, options):
    name = options.name
    uid = name
    species = []
    # Found to be the fastest way to strip non-alphanumeric characters
    # from a string in some post on StackOverflow
    pattern = re.compile(r'[\W]+')
    uid = pattern.sub('_', uid).strip('_')
    url = options.nexus
    with open(url_download(url, True), 'r') as fh:
        species = [line.strip(' (),').split(':')[0] for line in fh.readlines()]
    return name, uid.upper(), species


def get_maf_listing(maf_path):
    maf_files = []
    maf_url = urllib.parse.urlparse(maf_path)
    f = ftplib.FTP()
    f.connect(maf_url.netloc)
    f.login()
    listing = f.mlsd(maf_url.path)
    compressions = ['gz', 'bz2', 'zip']
    for name, facts in listing:
        skip = False
        if os.path.splitext(name)[-1].lstrip('.') not in compressions:
            skip = True
        if facts['type'] != 'file':
            skip = True
        for compression in compressions:
            for exclusion in ['_alt', '_random']:
                if name.endswith('%s.maf.%s' % (exclusion, compression)):
                    skip = True
                    break
        if name.startswith('chrUn'):
            skip = True
        if skip:
            continue
        maf_files.append(urllib.parse.urljoin(maf_path, name))
    f.close()
    return maf_files


def index_maf_files(maf_files, maf_path, options, params, target_directory):
    for maf_file in maf_files:
        maf_url = urllib.parse.urljoin(maf_path, maf_file)
        local_maf = url_download(maf_url, localpath=target_directory)
        index_command = ['maf_build_index.py', local_maf, local_maf + '.index']
        executor = subprocess.Popen(index_command)
        stdout, stderr = executor.communicate()


def main():
    parser = optparse.OptionParser()
    parser.add_option('-x', '--nexus', dest='nexus', action='store', type='string', help='URL for .nh')
    parser.add_option('-a', '--alignments', dest='alignments', action='store', type='string', help='URL for alignments')
    parser.add_option('-n', '--name', dest='name', action='store', type='string', help='Name')
    parser.add_option('-o', '--output', dest='output', action='store', type='string', help='Output')
    parser.add_option('-d', '--dbkey', dest='dbkey', action='store', type='string', help='dbkey')
    (options, args) = parser.parse_args()

    params = {}

    with open(options.output) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.makedirs(target_directory, exist_ok=True)

    display_name, uid, species_list = generate_metadata(params, options)
    maf_path = urllib.parse.urljoin(options.nexus, 'maf/')
    maf_files = get_maf_listing(maf_path)

    data_manager_entry = {
        'data_tables': {
            'indexed_maf_files': {
                'name': display_name,
                'dbkey': options.dbkey,  # This is needed for the output path
                'value': uid,
                'indexed_for': ','.join(species_list),
                'exists_in_maf': ','.join(species_list),
                'path': ','.join([maf_file.split('/')[-1] for maf_file in maf_files]),
            }
        }
    }

    # Fetch and index the MAFs
    index_maf_files(maf_files, maf_path, options, params, target_directory)
    with open(options.output, 'w') as fh:
        json.dump(data_manager_entry, fh, sort_keys=True)


if __name__ == "__main__":
    main()
