#!/usr/bin/env python

"""
Wrapper for syndiva.py
"""
import logging
import os
import shutil
import subprocess
import sys

import numpy

log = logging.getLogger(__name__)

assert sys.version_info[:2] >= (2, 4)


def stop_err(msg):
    sys.stderr.write("%s\n" % msg)
    sys.exit()


def __main__():
    # Parse Command Line

    fasta_file = sys.argv[1]
    pattern = sys.argv[2]
    restriction_site_5 = sys.argv[3]
    restriction_site_3 = sys.argv[4]
    extra_file_path = sys.argv[5] + "/"
    report = sys.argv[6]

    [ACDEFGHIKLMNPQRSTVWY\:][0-9][\-\*]
    authorized_pattern_letter = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                                 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', ':', '0',
                                 '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '*']
    authorized_restriction_site_letter = ['A', 'T', 'G', 'C']

    if numpy.array([letter in authorized_pattern_letter for letter in pattern]).sum() < len(pattern):
        sys.exit('1- Pattern Error running SynDivA: Pattern character does not match any allowed characters')
    if numpy.array([letter in authorized_restriction_site_letter for letter in restriction_site_5]).sum() < len(restriction_site_5):
        sys.exit('1- restriction site 5\' Error running SynDivA: restriction site 5\' character does not match any allowed characters')
    if numpy.array([letter in authorized_restriction_site_letter for letter in restriction_site_3]).sum() < len(restriction_site_3):
        sys.exit('1- restriction site 3\' Error running SynDivA: restriction site 3\' character does not match any allowed characters')
    try:  # for test - needs this don
        os.makedirs(extra_file_path)
    except Exception as e:
        sys.exit('1- Error running SynDivA ' + str(e))
    cmdline = 'python ./syndiva.py ' \
              '-i %s ' \
              '-o %s ' \
              '-p %s ' \
              '-5 %s ' \
              '-3 %s > /dev/null' % \
              (fasta_file,
               extra_file_path,
               pattern,
               restriction_site_5,
               restriction_site_3)
    try:
        proc = subprocess.Popen(args=cmdline,
                                shell=True,
                                stderr=subprocess.PIPE
                                )
        returncode = proc.wait()
        # get stderr, allowing for case where it's very large
        stderr = b''
        buffsize = 1048576
        try:
            while True:
                stderr += proc.stderr.read(buffsize)
                if not stderr or len(stderr) % buffsize != 0:
                    break
        except OverflowError:
            pass
        if returncode != 0:
            raise Exception(stderr)
    except Exception as e:
        stop_err('2 -Error running SynDivA ' + str(e))
    shutil.move(extra_file_path + "/syndiva_report.html", report)


if __name__ == "__main__":
    __main__()
