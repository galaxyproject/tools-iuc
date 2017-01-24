#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import re


def generate_index_html(dir_list, args):
    with open(args.html_file, 'w') as index_html_file:
        s = ""
        s += '<html>\n'
        s += '\t<head><title>PCoA beta diversity results</title></head>\n'
        s += '\t<body>\n'
        s += '\t\t<a href="http://www.qiime.org" target="_blank">'
        s += '<img src="http://qiime.org/_static/wordpressheader.png" '
        s += 'alt="www.qiime.org""/></a>\n'
        s += '\t\t<p>\n'
        s += '\t\t\tBeta diversity metrics\n'
        s += '\t\t\t<ul>\n'

        for directory in dir_list:
            regexp_result = re.search(
                '([a-zA-Z\_]*)_emperor_pcoa_plot',
                directory)
            metric = regexp_result.group(1)
            s += '\t\t\t\t<li>' + metric + ': '
            s += '<a href="' + directory
            s += '/index.html">PCoA results</a></td>\n'
            s += '\t\t\t\t</li>\n'

        s += '\t\t\t</ul>\n'
        s += '\t\t</p>\n'
        s += '\t</body>\n'
        s += '</html>\n'

        index_html_file.write(s)


def build_html(args):
    os.mkdir(args.html_dir)

    dir_list = [name for name in os.listdir(args.data_directory)
                if os.path.isdir(os.path.join(
                    args.data_directory,
                    name))]

    generate_index_html(dir_list, args)

    for directory in dir_list:
        input_path = os.path.join(args.data_directory, directory)
        cmd = 'cp -r ' + input_path + ' ' + args.html_dir
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_directory', required=True)
    parser.add_argument('--html_file', required=True)
    parser.add_argument('--html_dir', required=True)
    args = parser.parse_args()

    build_html(args)
