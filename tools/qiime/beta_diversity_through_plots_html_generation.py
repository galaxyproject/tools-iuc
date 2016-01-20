#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import re


def generate_index_html(dir_list, args):
    with open('index.html','w') as index_html_file:
        index_html_file.write('<html>\n')
        index_html_file.write('\t<head><title>PCoA beta diversity results</title></head>\n')
        index_html_file.write('\t<body>\n')
        index_html_file.write('\t\t<a href="http://www.qiime.org" target="_blank"><img src="http://qiime.org/_static/wordpressheader.png" alt="www.qiime.org""/></a>\n')
        index_html_file.write('\t\t<p>\n')
        index_html_file.write('\t\t\tBeta diversity metrics\n')
        index_html_file.write('\t\t\t<ul>\n')

        for directory in dir_list:
            regexp_result = re.match('([a-zA-Z]*)_emperor_pcoa_plot',directory)
            metric = regexp_result.group(1)
            index_html_file.write('\t\t\t\t<li>' + metric + ': ')
            index_html_file.write('<a href="' + directory + 
                '/index.html">PCoA results</a></td>\n')
            index_html_file.write('\t\t\t\t</li>\n')

        index_html_file.write('\t\t\t</ul>\n')
        index_html_file.write('\t\t</p>\n')
        index_html_file.write('\t</body>\n')
        index_html_file.write('</html>\n')

def build_html(args):
    os.mkdir(args.html_dir)

    dir_list = [ name for name in os.listdir(args.data_directory) if 
        os.path.isdir(os.path.join(args.data_directory, name)) ]

    generate_index_html(dir_list, args)

    os.system('cp -r index.html ' + args.html_dir)

    for directory in dir_list:
        os.system('cp -r ' + args.data_directory + '/' + directory + ' ' + 
            args.html_dir)
    
    os.system('mv ' + args.html_dir + '/index.html ' + args.html_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_directory', required=True)
    parser.add_argument('--html_file', required=True)
    parser.add_argument('--html_dir', required=True)
    args = parser.parse_args()

    build_html(args)
