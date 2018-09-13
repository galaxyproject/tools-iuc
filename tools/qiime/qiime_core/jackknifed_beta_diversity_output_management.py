#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os


def generate_index_html(dir_list, args):
    with open(args.html_file, 'w') as index_html_file:
        s = '<html>\n'
        s += '\t<head>\n'
        s += '\t\t<title>Jackknifed beta diversity results</title>\n'
        s += '\t</head>\n'
        s += '\t<body>\n'
        s += '\t\t<a href="http://www.qiime.org" target="_blank">'
        s += '<img src="http://qiime.org/_static/wordpressheader.png"'
        s += 'alt="www.qiime.org""/></a>\n'
        s += '\t\t<p>\n'
        s += '\t\t\tEmperor PCoA plot after jacknifed beta diversity\n'
        s += '\t\t\t<ul>\n'

        for directory in dir_list:
            if directory == 'rarefaction' or directory == 'unrarefied_bdiv':
                continue

            s += '\t\t\t\t<li><a href="' + directory
            s += '/index.html">PCoA plots</a></td></li>\n'

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

    os.mkdir(args.data_directory + '/pcoa')
    os.mkdir(args.data_directory + '/rare_dm')
    os.mkdir(args.data_directory + '/rare_upgma')
    os.mkdir(args.data_directory + '/upgma_cmp')
    os.mkdir(args.data_directory + '/rare_upgma_consensus')

    for directory in dir_list:
        if directory == 'rarefaction' or directory == 'unrarefied_bdiv':
            continue

        input_dir = os.path.join(args.data_directory, directory)

        input_path = os.path.join(input_dir, "emperor_pcoa_plots")
        output_path = os.path.join(args.html_dir, directory)
        cmd = 'cp -r ' + input_path + ' ' + output_path
        os.system(cmd)

        input_path = os.path.join(input_dir, "pcoa/*")
        output_path = os.path.join(args.data_directory, "pcoa")
        cmd = 'cp -r ' + input_path + ' ' + output_path
        os.system(cmd)

        input_path = os.path.join(input_dir, "rare_dm/*")
        output_path = os.path.join(args.data_directory, "rare_dm")
        cmd = 'cp -r ' + input_path + ' ' + output_path
        os.system(cmd)

        input_path = os.path.join(input_dir, "rare_upgma/*")
        output_path = os.path.join(args.data_directory, "rare_upgma")
        cmd = 'cp -r ' + input_path + ' ' + output_path
        os.system(cmd)

        input_path = os.path.join(input_dir, "upgma_cmp/*")
        output_path = os.path.join(args.data_directory, "upgma_cmp")
        cmd = 'cp -r ' + input_path + ' ' + output_path
        os.system(cmd)

        input_path = os.path.join(input_dir, "rare_upgma_consensus.tre")
        output_path = os.path.join(
            args.data_directory,
            "rare_upgma_consensus",
            directory + '.tre')
        cmd = 'cp ' + input_path + ' ' + output_path
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_directory', required=True)
    parser.add_argument('--html_file', required=True)
    parser.add_argument('--html_dir', required=True)
    args = parser.parse_args()

    build_html(args)
