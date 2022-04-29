#!/usr/bin/env python

import json
import os
import re
import sys
import tarfile
from urllib.request import urlretrieve


def main():
    # Read in given out_file and create target directory for file download
    with open(sys.argv[1]) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)

    # Process parameters for metadata and file download
    url = params['param_dict']['url'].rstrip("/") + "/" + params['param_dict']['file_name'].lstrip("/")
    m = re.search(r"(.*?)(merged|refseq)?_vep_(\d+?)_", params['param_dict']['file_name'])
    version = str(m.group(3))
    cache_type = m.group(2) if m.group(2) else "default"
    species = m.group(1).rstrip("_")
    display_name = f"{species.capitalize().replace('_', ' ')} {params['param_dict']['dbkey']} (V{version}{'' if cache_type == 'default' else ', ' + cache_type.capitalize()})"

    # Download and extract given cache archive, remove archive afterwards
    final_file, headers = urlretrieve(url, os.path.join(target_directory, params['param_dict']['file_name']))
    tar = tarfile.open(final_file, "r:gz")
    tar.extractall(target_directory)
    tar.close()
    os.remove(final_file)

    # Construct metadata for the new data table entry
    data_manager_dict = {
        'data_tables': {
            'vep_versioned_annotation_cache': [
                {
                    'value': params['param_dict']['file_name'].strip(".tar.gz"),
                    'dbkey': params['param_dict']['dbkey'],
                    'version': version,
                    'cachetype': cache_type,
                    'name': display_name,
                    'species': species,
                    'path': './%s' % params['param_dict']['file_name'].strip(".tar.gz")
                }
            ]
        }
    }

    # Save metadata to out_file
    with open(sys.argv[1], 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
