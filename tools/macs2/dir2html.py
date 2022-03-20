#!/usr/bin/env python

from __future__ import print_function

import os
import sys
from xml.sax.saxutils import escape


def make_table(directory):
    ret = ['<table class="fileList">\n']
    for filename in os.listdir(directory):
        ret.append('<tr><td class="file"><a href="%s">%s</a></td></tr>\n' % (filename, escape(filename).replace('MACS2_', '')))
    ret.append('</table>')
    return ''.join(ret)


def make_html(directory, stderr):
    return '\n'.join(['<html>'
                      '<head>',
                      '   <title>Additional output created by MACS2</title>',
                      '   <style type="text/css">',
                      '      table.fileList { text-align: left; }',
                      '      td.directory { font-weight: bold; }',
                      '      td.file { padding-left: 4em; }',
                      '   </style>',
                      '</head>',
                      '<body>',
                      '<h1>Additional Files:</h1>',
                      make_table(directory),
                      '<h3>Messages from MACS2:</h3>',
                      stderr.read().replace('\n', '<br>'),
                      '</body>',
                      '</html>'])


if __name__ == '__main__':
    if len(sys.argv) == 3:
        directory_path = sys.argv[1]
        stderr = open(sys.argv[2])
        print(make_html(directory_path, stderr))
    else:
        sys.exit('Two parameter expected: directory path and stderr path')
