"""
Wrapper for Prokka - Prokaryotic annotation tool
Author: Paolo Uva paolo dot uva at crs4 dot it
Date: February 14, 2013
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    op = optparse.OptionParser()
    op.add_option( '-f', '--fasta', dest="fasta", help="Fasta file with contigs")
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
    cl = ['prokka --outdir . --prefix prokka ', '%s' % (options.fasta)]
 
    # Run command
    dummy,tlog = tempfile.mkstemp(prefix='process_log')
    sout = open(tlog, 'w')
    p = subprocess.Popen(' '.join(cl), shell=True, stderr=sout, stdout=sout)
    retval = p.wait()
    sout.close()
    
    # Rename output files
    suffix = ['gff', 'gbk', 'fna', 'faa', 'ffn', 'sqn', 'fsa', 'tbl', 'err', 'log']
    try:
    	  for s in suffix:
    	  	shutil.move( 'prokka.' + s, getattr(options, s))
    except Exception, e:
    	  raise Exception, 'Error moving output file: ' + str( e )

if __name__=="__main__": __main__()
