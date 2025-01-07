#!/usr/bin/env python3

import argparse
import json
import sys
import tarfile

from io import BytesIO, StringIO
from pathlib import Path
from urllib.request import Request, urlopen

DATA_TABLE_NAME = 'clair3_models'


def find_latest_models():
    # based on the README.rst of the rerio repository as of 7 January 2025
    url = 'https://raw.githubusercontent.com/nanoporetech/rerio/refs/heads/master/README.rst'
    httprequest = Request(url)
    with urlopen(httprequest) as response:
        if response.status != 200:
            raise IOError(f'Failed to fetch the latest models: {response.status}')
        data = response.read().decode('utf-8')
        init_line_seen = False
        config_line_seen = False
        read_lines = False
        models = []
        break1 = 0
        for line in StringIO(data):
            if read_lines:
                if line.startswith('=========================='):
                    read_lines = False
                    break
                model = line[:break1 - 1]
                models.append(model)
            if config_line_seen and line.startswith('=========================='):
                break1 = line.find(' ')
                read_lines = True
                continue
            if init_line_seen and line.startswith('Config'):
                config_line_seen = True
                continue
            if line.startswith('Clair3 models for the following configurations are available:'):
                init_line_seen = True
                continue
        return models


def fetch_model(model_name):
    # the model files are tar gzipped, with a structure like:
    # model_name/pileup.index
    # model_name/full_alignment.index
    # and other files, with the key point being that the model_name becoomes the model_directory

    url = f'https://raw.githubusercontent.com/nanoporetech/rerio/refs/heads/master/clair3_models/{model_name}_model'
    httprequest = Request(url)
    with urlopen(httprequest) as response:
        if response.status != 200:
            raise IOError(f'Failed to fetch the model {model_name}: {response.status}')
        final_url = response.read().decode('utf-8').strip()
    httprequest = Request(final_url)
    with urlopen(httprequest) as response:
        if response.status != 200:
            raise IOError(f'Failed to fetch the model {model_name} from CDN URL {final_url}: {response.status}')
        data = response.read()
    return data


def unpack_model(data, outdir):
    with tarfile.open(fileobj=BytesIO(data), mode='r:*') as tar:
        tar.extractall(outdir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dm_filename', type=str, help='The filename of the data manager file to read parameters from and write outputs to')
    parser.add_argument('--download_latest', action='store_true', default=False, help='Download the latest models as per the rerio repository')
    parser.add_argument('--download_models', type=str, help='Comma separated list of models to download')
    args = parser.parse_args()

    # parameters to a data manager are passed in a JSON file (see https://docs.galaxyproject.org/en/latest/dev/data_managers.html) and
    # similarily a JSON file is created to pass the output back to Galaxy
    models = []
    if args.download_latest:
        models.extend(find_latest_models())
    if args.download_models:
        models.extend(args.download_models.split(','))

    if not models:
        sys.exit('No models to download, please specify either --download_latest or --download_models')

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)
    if 'extra_files_path' not in config.get('output_data', [{}])[0]:
        sys.exit('Please specify the output directory in the data manager configuration (the extra_files_path)')
    output_directory = config["output_data"]["extra_files_path"]
    if not Path(output_directory).exists():
        Path(output_directory).mkdir(parents=True)

    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = []

    for model in models:
        model_dir = Path(output_directory) / model
        # In the test below we assume that the contents of the model are uniquely identified by the model name. It is possible
        # that Oxford Nanopore will change the contents of the model tarball without changing the name. Hopefully this will never
        # happen.
        if model_dir.exists():
            print(f'Model {model} already exists, skipping', file=sys.stderr)
            continue
        data = fetch_model(model)
        unpack_model(data, output_directory)

        data_manager_dict["data_tables"][DATA_TABLE_NAME].append(
            dict(
                value=model,
                path=str(model_dir)
            )
        )

        with open(args.output_dm_filename, 'w') as fh:
            json.dump(data_manager_dict, fh, sort_keys=True, indent=4)
