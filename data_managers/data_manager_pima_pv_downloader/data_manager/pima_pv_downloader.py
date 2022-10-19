import argparse
import datetime
import json
import os
import sys
from urllib.request import Request, urlopen


def url_download(url, workdir):
    file_path = os.path.abspath(os.path.join(workdir, os.path.basename(url)))
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
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src:
            src.close()
    return file_path


def download(url, out_file):
    today = datetime.datetime.utcnow().strftime("%Y-%m-%d")

    with open(out_file) as fh:
        params = json.load(fh)

    workdir = params['output_data'][0]['extra_files_path']
    os.makedirs(workdir)
    file_path = url_download(url, workdir)
    name = '%s_%s' % (today, os.path.basename(file_path))

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = today
    data_manager_entry['name'] = name
    data_manager_entry['path'] = file_path
    data_manager_json["data_tables"]["pima_pv"] = data_manager_entry

    with open(out_file, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


parser = argparse.ArgumentParser()

parser.add_argument('--url', dest='url', help='URL to download plasmids_and_vectors.fasta file')
parser.add_argument('--out_file', dest='out_file', help='JSON output file')

args = parser.parse_args()

download(args.url, args.out_file)
