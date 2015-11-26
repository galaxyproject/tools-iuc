#!/usr/bin/env python

"""
Runs TreeVector on a newick file;
TODO: more documentation
"""

import optparse, os, shutil, subprocess, sys, tempfile, re, string

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='input', help='The sequence input file' )
    parser.add_option( '-s', '--shape', dest='shape', help='Branch shape' )
    parser.add_option( '-l', '--length', dest='length', help='Branch length' )
    parser.add_option( '-g', '--svg', dest='svg', help='Graph in SVG format' )
    parser.add_option( '-p', '--jarBin', dest='jarBin', default='', help='The path to where jars are stored' )
    parser.add_option( '-j', '--jarFile', dest='jarFile', help='The file name of the jar file to use')
    parser.add_option( '-x', '--jvmArgs', dest='jvmArgs', help='Java JVM arguments, e.g -Xmx250m')
    (options, args) = parser.parse_args()
    if options.jarBin == None:
       stop_err("Misssing option --jarBin")
    elif options.jarFile == None:
       stop_err("Misssing option --jarFile")
    elif options.input == None:
       stop_err("Misssing option --input")
    params = []
    props = []
    if options.jvmArgs != None:
        props.append(options.jvmArgs)
    if options.shape != None and options.shape != 'None':
       params.append('-%s' % options.shape)
    if options.length != None and options.length != 'None':
       params.append('-%s' % options.length)
    if options.svg != None and options.svg != 'None':
       params.append('-out %s' % options.svg)
    # make temp directory
    buffsize = 1048576
    tmp_dir = tempfile.mkdtemp()
    # print("tmp_dir %s" % tmp_dir)
    # generate commandline
    cmd = 'java %s -jar %s %s %s' % (' '.join(props), os.path.join( options.jarBin, options.jarFile ), options.input, ' '.join(params))
    # print >> sys.stderr, cmd 
    # need to nest try-except in try-finally to handle 2.4
    try:
        try:
            proc = subprocess.Popen( args=cmd, shell=True, stderr=subprocess.PIPE )
            returncode = proc.wait()
            stderr = proc.stderr.read()
            if returncode != 0:
                raise Exception, stderr
        except Exception, e:
            raise Exception, 'Error executing TeeVector. ' + str( e )
    except Exception, e:
        stop_err( 'TreeVector failed.\n' + str( e ) )

if __name__=="__main__": __main__()
