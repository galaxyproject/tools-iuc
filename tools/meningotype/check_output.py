#!/usr/bin/env python3

import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('meningotype_file', type=argparse.FileType())
    args = parser.parse_args()
    for line in args.meningotype_file:
        if 'ERROR:' in line:
            sys.exit('Found error in meningotype output: ' + line.strip())
