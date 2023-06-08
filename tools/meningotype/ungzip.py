#!/usr/bin/env python3

import argparse
import gzip

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_input')
    parser.add_argument('dataset_output', type=argparse.FileType('wb'))
    args = parser.parse_args()
    
    blocksize = 1024*1000*100;  # 100 MB buffer
    input_file = gzip.open(args.dataset_input, "rb")
    data = input_file.read(blocksize)
    while len(data) != 0:
        args.dataset_output.write(data)
        data = input_file.read(blocksize)
