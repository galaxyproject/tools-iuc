import argparse
import datetime
import json
import os
import sys
import uuid

import requests
from requests.exceptions import ContentDecodingError


def url_download(url):
    """Attempt to download gene annotation file from a given url
    :param url: full url to gene annotation file
    :type url: str.
    :returns: name of downloaded gene annotation file
    :raises: ContentDecodingError, IOError
    """
    response = requests.get(url=url, stream=True)

    # Generate file_name
    file_name = response.url.split("/")[-1]

    block_size = 10 * 1024 * 1024  # 10MB chunked download
    with open(file_name, 'w+') as f:
        try:
            # Good to note here that requests' iter_content() will
            # automatically handle decoding "gzip" and "deflate" encoding
            # formats
            for buf in response.iter_content(block_size):
                f.write(buf)
        except (ContentDecodingError, IOError) as e:
            sys.stderr.write("Error occured downloading reference file: %s"
                             % e)
            os.remove(file_name)

    return file_name


def main():
    parser = argparse.ArgumentParser(description='Create data manager JSON.')
    parser.add_argument('--out', dest='output', action='store',
                        help='JSON filename')
    parser.add_argument('--name', dest='name', action='store',
                        default=uuid.uuid4(), help='Data table entry unique ID'
                        )
    parser.add_argument('--url', dest='url', action='store',
                        help='Url to download gtf file from')

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

    with open(os.path.join(args.output), "w+") as f:
        f.write(json.dumps(data_manager_entry))


if __name__ == '__main__':
    main()
