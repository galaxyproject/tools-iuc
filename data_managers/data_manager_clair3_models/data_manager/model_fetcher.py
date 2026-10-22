#!/usr/bin/env python3

import argparse
import json
import sys
import tarfile
from hashlib import sha256
from io import BytesIO, StringIO
from pathlib import Path
from urllib.error import HTTPError
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
        latest_seen = False
        config_line_seen = False
        read_lines = False
        models = []
        # the file that we are parsing has a section that looks like this:
        # Clair3 Models
        # -------------

        # Clair3 models for the following configurations are available:

        # Latest:

        # ========================== =================== =======================
        # Config                     Chemistry           Dorado basecaller model
        # ========================== =================== =======================
        # r1041_e82_400bps_sup_v500  R10.4.1 E8.2 (5kHz) v5.0.0 SUP
        # r1041_e82_400bps_hac_v500  R10.4.1 E8.2 (5kHz) v5.0.0 HAC
        # r1041_e82_400bps_sup_v410  R10.4.1 E8.2 (4kHz) v4.1.0 SUP
        # r1041_e82_400bps_hac_v410  R10.4.1 E8.2 (4kHz) v4.1.0 HAC
        # ========================== =================== =======================
        #
        # and the aim is to extract the list of model names from the table by successfully looking for
        # "Clair3 Models", then "Latest:", then "Config" and then "=====" and then reading the lines until
        # the next "=====" is encountered
        for line in StringIO(data):
            if read_lines:
                if line.startswith('====='):
                    read_lines = False
                    break
                model = line.split()[0]
                models.append(model)
            if config_line_seen and line.startswith('====='):
                read_lines = True
                continue
            if init_line_seen and line.startswith('Latest:'):
                latest_seen = True
                continue
            if latest_seen and line.startswith('Config'):
                config_line_seen = True
                continue
            if line.startswith('Clair3 Models'):
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
    try:
        # urlopen throws a HTTPError if it gets a 404 status (and perhaps other non-200 status?)
        with urlopen(httprequest) as response:
            if response.status != 200:
                raise IOError(f'Failed to fetch the model {model_name}: {response.status}')
            final_url = response.read().decode('utf-8').strip()
        httprequest = Request(final_url)
    except HTTPError as e:
        raise IOError(f'Failed to fetch the model {model_name}: {e}')

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
    parser.add_argument('--known_models', type=str, help='List of models already known in the Galaxy data table')
    parser.add_argument('--sha256_sums', type=str, help='List of sha256sums of the models already known in the Galaxy data table')
    parser.add_argument('--download_latest', action='store_true', default=False, help='Download the latest models as per the Rerio repository')
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

    with open(args.dm_filename) as fh:
        config = json.load(fh)
    if 'extra_files_path' not in config.get('output_data', [{}])[0]:
        sys.exit('Please specify the output directory in the data manager configuration (the extra_files_path)')
    output_directory = config["output_data"][0]["extra_files_path"]
    if not Path(output_directory).exists():
        Path(output_directory).mkdir(parents=True)

    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][DATA_TABLE_NAME] = []

    known_models = set(args.known_models.split(',')) if args.known_models else set()
    model_to_sha256 = {}
    if args.known_models:
        sha256_sums = args.sha256_sums.split(',')
        for (i, model) in enumerate(known_models):
            model_to_sha256[model] = sha256_sums[i]

    for model in models:
        model_dir = Path(output_directory) / model
        # The data table cannot handle duplicate entries, so we skip models that are already in the data table
        if model in known_models:
            print(f'Model {model} already exists, skipping', file=sys.stderr)
            continue
        data = fetch_model(model)
        sha256sum = sha256(data).hexdigest()

        # Since we skip models that are already known we cannot test the sha256sum here. This code is retained to illustrate that an
        # alternative logic would be to download the model each time and check if the sha256sum matches what is already known. Hopefully
        # ONT does not update the models while keeping the same name, so this is not needed. The sha256sum is stored in the data table
        # in case it is needed in the future.
        # if model in model_to_sha256 and sha256sum != model_to_sha256[model]:
        #    sys.exit(f'Model {model} already exists with a different sha256sum {model_to_sha256[model]}. This is a serious error, inform the Galaxy admin')

        unpack_model(data, output_directory)

        data_manager_dict["data_tables"][DATA_TABLE_NAME].append(
            dict(
                value=model,
                platform="ont",
                sha256=sha256sum,
                path=str(model_dir),
                source="rerio"
            )
        )

    with open(args.dm_filename, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True, indent=4)
