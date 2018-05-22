#!/usr/bin/env python

# This file is derived from code in the Galaxy project that, but due to historical
# reasons reflecting time developed outside of the Galaxy Project, this file is under
# the MIT license.
#
# The MIT License (MIT)
# Copyright (c) 2012,2013,2014,2015,2016 Peter Cock
# Copyright (c) 2012 Edward Kirton
# Copyright (c) 2013 Nicola Soranzo
# Copyright (c) 2014 Bjoern Gruening
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE

import re
import shutil
import sys
xmlfiles = sys.argv[2:]
xmlout = sys.argv[1]


def merge(split_files, output_file):
    iter_num_re = re.compile(' *<Iteration_iter-num>\d+</Iteration_iter-num>')
    query_ID_re = re.compile(' *<Iteration_query-ID>Query_\d+</Iteration_query-ID>')
    if len(split_files) == 1:
        shutil.copy(split_files[0], output_file)
        # TODO return merge(split_files, output_file)
    out = open(output_file, "w")
    h = None
    iter_num = 2  # we only start using this from the 2nd file onwards
    for f in split_files:
        h = open(f)
        # body = False
        first_line = False
        header = h.readline()
        if not header:
            out.close()
            h.close()
            raise ValueError("BLAST XML file %s was empty" % f)
        if header.strip() != '<?xml version="1.0"?>':
            out.write(header)  # for diagnosis
            out.close()
            h.close()
            raise ValueError("%s is not an XML file!" % f)
        line = h.readline()
        header += line
        # check that for BLAST doctype
        if line.strip() not in [
                '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
                '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">']:
            out.write(header)  # for diagnosis
            out.close()
            h.close()
            raise ValueError("%s is not a BLAST XML file!" % f)
        # read in header (part before the first <Iteration>)
        while True:
            line = h.readline()
            if not line:
                out.write(header)  # for diagnosis
                out.close()
                h.close()
                raise ValueError("BLAST XML file %s ended prematurely" % f)
            header += line
            if "<Iteration>" in line:
                break
            if len(header) > 10000:
                # Something has gone wrong, don't load too much into memory!
                # Write what we have to the merged file for diagnostics
                out.write(header)
                out.close()
                h.close()
                raise ValueError("BLAST XML file %s has too long a header!" % f)
        if "<BlastOutput>" not in header:
            out.close()
            h.close()
            raise ValueError("%s is not a BLAST XML file:\n%s\n..." % (f, header))
        if f == split_files[0]:
            out.write(header)
            old_header = header
        elif old_header[:300] != header[:300]:
            # Enough to check <BlastOutput_program> and <BlastOutput_version> match
            out.close()
            h.close()
            raise ValueError("BLAST XML headers don't match for %s and %s - have:\n%s\n...\n\nAnd:\n%s\n...\n" %
                             (split_files[0], f, old_header[:300], header[:300]))
        else:
            first_line = True
        for line in h:
            if first_line or "<Iteration>" in line:
                if "<Iteration>" in line:
                    line = next(h)
                out.write("<Iteration>\n")
                first_line = False
                iter_num_match = iter_num_re.match(line)
                if iter_num_match is not None:
                    out.write('  <Iteration_iter-num>{}</Iteration_iter-num>\n'.format(iter_num))
                else:
                    out.write(line)
                line = next(h)
                query_ID_match = query_ID_re.match(line)
                if query_ID_match is not None:
                    out.write('  <Iteration_query-ID>Query_{}</Iteration_query-ID>\n'.format(iter_num))
                else:
                    out.write(line)
                iter_num += 1
                continue

            if "</BlastOutput_iterations>" in line:
                break
            out.write(line)
        h.close()
    out.write("  </BlastOutput_iterations>\n")
    out.write("</BlastOutput>\n")
    out.close()


merge(xmlfiles, xmlout)
