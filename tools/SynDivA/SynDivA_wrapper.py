#!/usr/bin/env python

"""
Wrapper for SynDivA.py
"""
import pkg_resources
import logging, os, string, sys, tempfile, glob, shutil, types, urllib
import shlex, subprocess
from optparse import OptionParser, OptionGroup
from stat import *


log = logging.getLogger( __name__ )

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write("%s\n" % msg)
    sys.exit()

def __main__():
    #Parse Command Line
    s = 'SynDivA_wrapper.py:  argv = %s\n' % (sys.argv)
    argcnt = len(sys.argv)
    fasta_file = sys.argv[1]
    pattern = sys.argv[2]
    restriction_site_5 = sys.argv[3]
    restriction_site_3 = sys.argv[4]
    install_dir = sys.argv[5]
    extra_file_path = sys.argv[6]+"/"
    report = sys.argv[7]
    tmp_file_path = sys.argv[8]
    tool_file_path = sys.argv[9]+"/"
    try:# for test - needs this done
        os.makedirs(extra_file_path)
    except Exception as e:
        stop_err('1- Error running SynDivA ' + str(e))
    cmdline = 'python %sSynDivA.py -i %s -o %s -p %s -5 %s -3 %s > /dev/null' % (tool_file_path, fasta_file, extra_file_path, pattern, restriction_site_5, restriction_site_3)
    try:
        proc = subprocess.Popen(args=cmdline, shell=True, stderr=subprocess.PIPE)
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
    png_path = os.path.join(extra_file_path,'distrib.png')
    shutil.move(extra_file_path+"/SynDivA_report.html", report)
    #rval = ['<html><head><title>SynDivA Galaxy Composite Dataset </title></head><p/>']
    #rval.append('<div>%s<p/></div>' % (cmdline) )
    #rval.append('<div>This composite dataset is composed of the following files:<p/><ul>')
    #rval.append( '<li><a href="%s" type="text/plain">%s </a>%s</li>' % (png_path,'Sequences','Sequences' ) )
    #rval.append( '</ul></div></html>' )
    #f = file(html_file,'w')
    #f.write("\n".join( rval ))
    #f.write('\n')
    #f.close()

if __name__ == "__main__": __main__()
