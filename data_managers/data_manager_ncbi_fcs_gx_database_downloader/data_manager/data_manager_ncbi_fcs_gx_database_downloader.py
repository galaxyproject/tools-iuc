#!/usr/bin/env python

import argparse
import json
import os
import sys
import subprocess



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--manifest', required=True)
    parser.add_argument('--local_cache_dir', required=True)
    parser.add_argument('--tool_cache_dir', required=True)

    args = parser.parse_args()

    print (args)


# ncbi_fcs_gx_databases.loc
# 
# 
# 
# caching | running ; tag ;  src manifest (url or path) ; dst dir
# 
# 
# 
# 
# 
# 
# local_dir 
# 
# 
# 
# 
# manager
# 
# 
#           "context": "",
#           "tag": "",
#           "src_manifest_dir": "",
#           "dst_dir": ""
