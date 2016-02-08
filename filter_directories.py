#!/usr/bin/env python

import os
import sys

with open('.tt_blacklist') as handle:
    bl = [tool.strip() for tool in handle]

for directory in sys.stdin:
    directory = directory.strip()
    if directory.startswith('packages') or directory.startswith('data_managers'):
        continue
    while directory:
        if os.path.exists( os.path.join(directory, '.shed.yml') ) or os.path.exists( os.path.join(directory, '.shed.yaml') ):
            if directory not in bl:
                print(directory)
            break
        else:
            directory = os.path.dirname( directory )
