#!/usr/bin/env python

import json
import os
import tarfile

# rely on the fact that pangolin itself uses the requests module
import requests

response = requests.get("https://api.github.com/repos/cov-lineages/pangoLEARN/releases/latest")
if response.status_code == 200:
    details = json.loads(response.text)
    response = requests.get(details["tarball_url"])
    if response.status_code == 200:
        open("pangolearn.tgz", "wb").write(response.content)
        tf = tarfile.open("pangolearn.tgz")
        pl_path = tf.next().name
        tf.extractall()
        os.rename(pl_path + '/pangoLEARN', 'datadir')
    else:
        response.raise_for_status()
else:
    response.raise_for_status()


