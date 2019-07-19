import argparse
import datetime
import json
import os
import shutil
import tarfile
import zipfile
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen


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


def main(args):
    workdir = os.path.join(os.getcwd(), 'taxonomy')
    url_download(args.url, workdir)
    data_manager_entry = {}
    data_manager_entry['value'] = args.name.lower()
    data_manager_entry['name'] = args.name
    data_manager_entry['path'] = '.'
    data_manager_json = dict(data_tables=dict(ncbi_taxonomy=data_manager_entry))
    params = json.loads(open(args.output).read())
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    output_path = os.path.abspath(os.path.join(os.getcwd(), 'taxonomy'))
    for filename in os.listdir(workdir):
        shutil.move(os.path.join(output_path, filename), target_directory)
    open(args.output, 'w').write(json.dumps(data_manager_json))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--out', dest='output', action='store', help='JSON filename')
    parser.add_argument('--name', dest='name', action='store', default=str(datetime.date.today()), help='Data table entry unique ID')
    parser.add_argument('--url', dest='url', action='store', default='ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz', help='Download URL')
    args = parser.parse_args()

    main(args)
