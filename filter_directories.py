#!/usr/bin/env python

import os
import sys

for directory in sys.stdin:
    directory = directory.strip()
    while directory:
        if os.path.exists( os.path.join(directory, '.shed.yml') ) or os.path.exists( os.path.join(directory, '.shed.yaml') ):
            print(directory)
            break
        else:
            directory = os.path.dirname( directory )
