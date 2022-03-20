import argparse
import datetime
import json
import os
import shutil
import tarfile
import zipfile
from urllib.request import Request, urlopen


def url_download(url, workdir):
    file_path = os.path.join(workdir, 'download.dat')
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(file_path, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    finally:
        if src:
            src.close()
    if tarfile.is_tarfile(file_path):
        fh = tarfile.open(file_path, 'r:*')
    elif zipfile.is_zipfile(file_path):
        fh = zipfile.ZipFile(file_path, 'r')
    else:
        return
    fh.extractall(workdir)
    os.remove(file_path)


def download_name_maps(url, workdir, partial):

    if partial:
        map_files = [
            'pdb.accession2taxid.gz',
        ]
    else:
        map_files = [
            'dead_nucl.accession2taxid.gz',
            'dead_prot.accession2taxid.gz',
            'dead_wgs.accession2taxid.gz',
            'nucl_gb.accession2taxid.gz',
            'nucl_wgs.accession2taxid.gz',
            'pdb.accession2taxid.gz',
            'prot.accession2taxid.gz',
            'prot.accession2taxid.FULL.gz'
        ]

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    for map in map_files:
        src = "{}{}".format(url, map)
        dest = os.path.join(workdir, map)

        print("Downloading taxonomy accession2taxid file from {} to {}".format(src, dest))

        try:
            req = Request(src)
            src = urlopen(req)
            with open(dest, 'wb') as dst:
                while True:
                    chunk = src.read(2**10)
                    if chunk:
                        dst.write(chunk)
                    else:
                        break
        finally:
            if src:
                src.close()


def move_files_to_final_dir(workdir, target_directory, copy=False):
    for filename in os.listdir(workdir):
        if copy:
            shutil.copy(os.path.join(workdir, filename), target_directory)
        else:
            shutil.move(os.path.join(workdir, filename), target_directory)


def main(args):
    workdir = os.path.abspath(os.path.join(os.getcwd(), 'taxonomy'))
    url_download(args.url, workdir)

    data_manager_entry = {}
    data_manager_entry['value'] = args.name.lower()
    data_manager_entry['name'] = args.name
    data_manager_entry['path'] = '.'
    data_manager_json = dict(data_tables=dict(ncbi_taxonomy=data_manager_entry))

    with open(args.output) as fh:
        params = json.load(fh)

    if args.name_maps:
        workdir_a2t = os.path.join(os.getcwd(), 'accession2taxid')
        download_name_maps("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/", workdir_a2t, args.partial)

        target_directory_a2t = os.path.join(params['output_data'][0]['extra_files_path'], "accession2taxid")
        os.makedirs(target_directory_a2t)
        move_files_to_final_dir(workdir_a2t, target_directory_a2t)

        # Also copy taxonomy data to accession2taxid dir
        move_files_to_final_dir(workdir, target_directory_a2t, copy=True)

        data_manager_json['data_tables']['ncbi_accession2taxid'] = data_manager_entry

    target_directory_tax = os.path.join(params['output_data'][0]['extra_files_path'], "taxonomy")
    os.makedirs(target_directory_tax)

    move_files_to_final_dir(workdir, target_directory_tax)

    with open(args.output, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--out', dest='output', action='store', help='JSON filename')
    parser.add_argument('--name', dest='name', action='store', default=str(datetime.date.today()), help='Data table entry unique ID')
    parser.add_argument('--url', dest='url', action='store', default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz', help='Download URL')
    parser.add_argument('--name-maps', dest='name_maps', action='store_true', help='')
    parser.add_argument('--partial', dest='partial', action='store_true', help='Only download a small subset of data (for testing)')
    args = parser.parse_args()

    main(args)
