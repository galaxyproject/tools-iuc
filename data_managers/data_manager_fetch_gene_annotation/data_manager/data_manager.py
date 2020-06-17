# -*- coding: utf-8 -*-
import argparse
import bz2
import datetime
import gzip
import json
import os
import shutil
import sys
import uuid
import zipfile


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
        return data.startswith(self.magic)

    def open(self):
        return None


class ZIPFile(CompressedFile):
    magic = '\x50\x4b\x03\x04'
    file_type = 'zip'
    mime_type = 'compressed/zip'

    def open(self):
        return zipfile.ZipFile(self.f)


class BZ2File(CompressedFile):
    magic = '\x42\x5a\x68'
    file_type = 'bz2'
    mime_type = 'compressed/bz2'

    def open(self):
        return bz2.BZ2File(self.f)


class GZFile(CompressedFile):
    magic = '\x1f\x8b\x08'
    file_type = 'gz'
    mime_type = 'compressed/gz'

    def open(self):
        return gzip.GzipFile(self.f)


# factory function to create a suitable instance for accessing files
def get_compressed_file(filename):
    with open(filename, 'rb') as f:
        start_of_file = f.read(1024)
        f.seek(0)
        for cls in (ZIPFile, BZ2File, GZFile):
            if cls.is_magic(start_of_file):
                f.close()
                return cls(filename)

        return None


try:
    # For Python 3.0 and later
    from urllib.request import urlretrieve
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib import urlretrieve


def url_download(url):
    """Attempt to download gene annotation file from a given url
    :param url: full url to gene annotation file
    :type url: str.
    :returns: name of downloaded gene annotation file
    :raises: ContentDecodingError, IOError
    """

    # Generate file_name
    file_name = url.split('/')[-1]

    try:
        # download URL (FTP and HTTP work, probably local and data too)
        urlretrieve(url, file_name)

        # uncompress file if needed
        cf = get_compressed_file(file_name)
        if cf is not None:
            uncompressed_file_name = os.path.splitext(file_name)[0]
            with open(uncompressed_file_name, 'w+') as uncompressed_file:
                shutil.copyfileobj(cf.accessor, uncompressed_file)
            os.remove(file_name)
            file_name = uncompressed_file_name
    except IOError as e:
        sys.stderr.write('Error occured downloading reference file: %s' % e)
        os.remove(file_name)
    return file_name


def main():
    parser = argparse.ArgumentParser(description='Create data manager JSON.')
    parser.add_argument('--out', dest='output', action='store', help='JSON filename')
    parser.add_argument('--name', dest='name', action='store', default=uuid.uuid4(), help='Data table entry unique ID')
    parser.add_argument('--url', dest='url', action='store', help='Url to download gtf file from')

    args = parser.parse_args()

    work_dir = os.getcwd()

    # Attempt to download gene annotation file from given url
    gene_annotation_file_name = url_download(args.url)

    # Update Data Manager JSON and write to file
    data_manager_entry = {
        'data_tables': {
            'gff_gene_annotations': {
                'value': str(datetime.datetime.now()),
                'dbkey': str(args.name),
                'name': gene_annotation_file_name,
                'path': os.path.join(work_dir, gene_annotation_file_name)
            }
        }
    }

    with open(os.path.join(args.output), 'w+') as f:
        f.write(json.dumps(data_manager_entry))


if __name__ == '__main__':
    main()
