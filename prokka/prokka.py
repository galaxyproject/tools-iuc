# -*- coding: utf-8 -*-
"""
Wrapper for Prokka - Prokaryotic annotation tool
Author: Paolo Uva paolo dot uva at crs4 dot it
Date: February 14, 2013
Update: March 14, 2013 - Added more options
"""

import optparse
import shutil
import subprocess
import sys

def stop_err( msg ):
    """ Print error message and exit """
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    op = optparse.OptionParser()
    op.add_option( '', '--fasta',     dest="fasta",                         help="Fasta file with contigs")
    op.add_option( '', '--kingdom',   dest="kingdom",   default="Bacteria", help="Kingdom")
    op.add_option( '', '--gram',      dest="gram",      default="None",     help="Gram")
    op.add_option( '', '--minContig', dest="minContig", default="200",      help="Minimun contig size")
    op.add_option( '', '--rfam',      dest="rfam",      default="false",    help="Enable searching for ncRNAs")
    op.add_option( '', '--centre',    dest="centre",    default="CRS4",     help="Sequencing centre")
    op.add_option( '', '--gff', dest="gff", help="This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV")
    op.add_option( '', '--gbk', dest="gbk", help="This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence")
    op.add_option( '', '--fna', dest="fna", help="Nucleotide FASTA file of the input contig sequences")
    op.add_option( '', '--faa', dest="faa", help="Protein FASTA file of the translated CDS sequences")
    op.add_option( '', '--ffn', dest="ffn", help="Nucleotide FASTA file of all the annotated sequences, not just CDS")
    op.add_option( '', '--sqn', dest="sqn", help="An ASN1 format Sequin file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc")
    op.add_option( '', '--fsa', dest="fsa", help="Nucleotide FASTA file of the input contig sequences, used by tbl2asn to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines")
    op.add_option( '', '--tbl', dest="tbl", help="Feature Table file, used by tbl2asn to create the .sqn file")
    op.add_option( '', '--err', dest="err", help="Unacceptable annotations - the NCBI discrepancy report")
    op.add_option( '', '--log', dest="log", help="Contains all the output that Prokka produced during its run")

    (options, args) = op.parse_args()

    # Build command
    if options.gram == 'None' :
        gram_flag = ''
    else:
        gram_flag = '--gram %s' % (options.gram)
    if options.rfam == 'false' :
        rfam_flag = ''
    else:
        rfam_flag = '--rfam'
    
    cl = 'prokka --force --outdir . --prefix prokka --kingdom %s --minContig %s --centre %s %s %s %s' % (options.kingdom, options.minContig, options.centre, gram_flag, rfam_flag, options.fasta)
    print '\nProkka command to be executed: \n %s' % ( cl )

    # Run command
    try:
        if options.log:
            sout = open(options.log, 'w')
        else:
            sout = sys.stdout
        retcode = subprocess.call(cl, stdout=sout, stderr=subprocess.STDOUT, shell=True) # need to redirect stderr because prokka writes many logging info there
        if sout != sys.stdout:
            sout.close()
        
        if retcode != 0:
            stop_err("Execution of Prokka terminated with return code: %i" % retcode)
    except Exception, e:
        stop_err('Execution of Prokka failed:: ' + str(e))
    
    # Rename output files
    suffix = ['gbk', 'fna', 'faa', 'ffn', 'sqn', 'fsa', 'tbl', 'err', 'gff']
    try:
        for s in suffix:
            shutil.move( 'prokka.' + s, getattr(options, s))
    except Exception, e:
        stop_err('Error moving output file: ' + str(e))

if __name__ == "__main__":
    __main__()
