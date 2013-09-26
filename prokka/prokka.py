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


def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('--cpus', dest='cpus', type='int', help='Number of CPUs to use [0=all]')
    parser.add_option('--fasta', dest='fasta', help='FASTA file with contigs')
    parser.add_option('--kingdom', dest='kingdom', choices=['Archaea', 'Bacteria', 'Viruses'], default='Bacteria', help='Kingdom')
    parser.add_option('--mincontig', dest='mincontig', type='int', help='Minimun contig size')
    parser.add_option('--rfam', action="store_true", dest="rfam", help="Enable searching for ncRNAs")
    parser.add_option('--centre', dest="centre", default="CRS4", help="Sequencing centre")
    parser.add_option('--gff', dest="gff", help="This is the master annotation in GFF3 format, containing both sequences and annotations")
    parser.add_option('--gbk', dest="gbk", help="This is a standard GenBank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-GenBank, with one record for each sequence")
    parser.add_option('--fna', dest="fna", help="Nucleotide FASTA file of the input contig sequences")
    parser.add_option('--faa', dest="faa", help="Protein FASTA file of the translated CDS sequences")
    parser.add_option('--ffn', dest="ffn", help="Nucleotide FASTA file of all the annotated sequences, not just CDS")
    parser.add_option('--sqn', dest="sqn", help="An ASN1 format Sequin file for submission to GenBank. It needs to be edited to set the correct taxonomy, authors, related publication, etc.")
    parser.add_option('--fsa', dest="fsa", help="Nucleotide FASTA file of the input contig sequences, used by tbl2asn to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines")
    parser.add_option('--tbl', dest="tbl", help="Feature Table file, used by tbl2asn to create the .sqn file")
    parser.add_option('--err', dest="err", help="Unacceptable annotations - the NCBI discrepancy report")
    parser.add_option('--txt', dest='txt', help='Statistics relating to the annotated features found')
    parser.add_option('--log', dest="log", help="Contains all the output that Prokka produced during its run")
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')

    # Build command
    cpus = "--cpus %d" % (options.cpus) if options.cpus is not None else ''
    rfam = '--rfam' if options.rfam else ''
    mincontig = "--mincontig %d" % options.mincontig if options.mincontig is not None else ''

    cl = "prokka --force --outdir . --prefix prokka --kingdom %s %s --centre %s %s %s %s" % (options.kingdom, mincontig, options.centre, rfam, cpus, options.fasta)
    print '\nProkka command to be executed:\n %s' % cl

    # Run command
    log = open(options.log, 'w') if options.log else sys.stdout
    try:
        subprocess.check_call(cl, stdout=log, stderr=subprocess.STDOUT, shell=True) # need to redirect stderr because prokka writes many logging info there
    finally:
        if log != sys.stdout:
            log.close()

    # Rename output files
    suffix = ['gff', 'gbk', 'fna', 'faa', 'ffn', 'sqn', 'fsa', 'tbl', 'err', 'txt']
    for s in suffix:
        shutil.move('prokka.' + s, getattr(options, s))

if __name__ == "__main__":
    __main__()
