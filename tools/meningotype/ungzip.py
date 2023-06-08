#!/usr/bin/env python3

import argparse
import gzip

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_input')
    parser.add_argument('dataset_name')
    args = parser.parse_args()
    
    blocksize = 1024*1000*100;  # 100 MB buffer
    input_file = gzip.open(args.dataset_input, "rb")
    output_filename = args.dataset_name.replace('.gz', '')
    output_file = open(output_filename, "wb")

    data = input_file.read(blocksize)
    while len(data) != 0:
        output_file.write(data)
        data = input_file.read(blocksize)
