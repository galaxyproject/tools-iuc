#!/usr/bin/env python

import argparse
import json
import os
import sys
import tarfile

# rely on the fact that pangolin itself uses the requests module
import requests
from requests.adapters import HTTPAdapter

parser = argparse.ArgumentParser(description='download latest pangoLEARN database')
parser.add_argument('--timeout', type=float, default=60.0, help="Timeout for connecting or downloading the pangoLEARN database (in seconds)")
parser.add_argument('--max_retries', type=int, default=5, help="How many times to retry fetching the pangoLEARN database")
args = parser.parse_args()

try:
    s = requests.Session()
    s.mount('https://', HTTPAdapter(max_retries=args.max_retries))
    response = s.get(
        "https://api.github.com/repos/cov-lineages/pangoLEARN/releases/latest"
    )
    if response.status_code == 200:
        details = json.loads(response.text)
        response = s.get(details["tarball_url"], timeout=args.timeout)
        if response.status_code == 200:
            with open("pangolearn.tgz", "wb") as handle:
                handle.write(response.content)
            tf = tarfile.open("pangolearn.tgz")
            pl_path = tf.next().name
            tf.extractall()
            tf.close()
            os.rename(os.path.join(pl_path, "pangoLEARN"), "datadir")
        else:
            response.raise_for_status()
    else:
        response.raise_for_status()
except (requests.ConnectionError, requests.HTTPError, requests.TooManyRedirects) as e:
    sys.exit('Failed to download pangoLEARN database: {}'.format(e))
except requests.Timeout:
    sys.exit('Timeout while trying to download pangoLEARN database')
