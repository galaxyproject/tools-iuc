#purpose: macs2 python wrapper
#author: Ziru Zhou
#date: November, 2012

import os
import sys
import subprocess
import tempfile
import shutil
import glob
import gzip
from galaxy import eggs
import pkg_resources
pkg_resources.require( "simplejson" )
import simplejson

CHUNK_SIZE = 1024

#==========================================================================================
#functions
#==========================================================================================
def gunzip_cat_glob_path( glob_path, target_filename, delete = False ):
    out = open( target_filename, 'wb' )
    for filename in glob.glob( glob_path ):
        fh = gzip.open( filename, 'rb' )
        while True:
            data = fh.read( CHUNK_SIZE )
            if data:
                out.write( data )
            else:
                break
        fh.close()
        if delete:
            os.unlink( filename )
    out.close()

def xls_to_interval( xls_file, interval_file, header = None ):
    out = open( interval_file, 'wb' )
    if header:
        out.write( '#%s\n' % header )
    wrote_header = False
    #From macs readme: Coordinates in XLS is 1-based which is different with BED format.
    for line in open( xls_file ):
        #keep all existing comment lines
        if line.startswith( '#' ):
            out.write( line )
        #added for macs2 since there is an extra newline 
        elif line.startswith( '\n' ):
            out.write( line )
        elif not wrote_header:
            out.write( '#%s' % line )
        #print line
            wrote_header = True
        else:
            fields = line.split( '\t' )
            if len( fields ) > 1:
                fields[1] = str( int( fields[1] ) - 1 )
            out.write( '\t'.join( fields ) )
    out.close()

#==========================================================================================
#main
#==========================================================================================
def main():
    #take in options file and output file names
    options = simplejson.load( open( sys.argv[1] ) )
    outputs = simplejson.load( open( sys.argv[2] ) )

    #=================================================================================
    #parse options and execute macs2
    #=================================================================================
    #default inputs that are in every major command
    experiment_name = '_'.join( options['experiment_name'].split() ) #save experiment name here, it will be used by macs for some file names
    cmdline = "macs2 %s -t %s" % ( options['command'], ",".join( options['input_chipseq'] ) )
    if options['input_control']:
        cmdline = "%s -c %s" % ( cmdline, ",".join( options['input_control'] ) )

    #=================================================================================
    if options['command'] == "callpeak":
        output_bed = outputs['output_bed_file']
        output_extra_html = outputs['output_extra_file']
        output_extra_path = outputs['output_extra_file_path']
        output_peaks =  outputs['output_peaks_file']
        output_narrowpeaks = outputs['output_narrowpeaks_file']
        output_xls_to_interval_peaks_file = outputs['output_xls_to_interval_peaks_file']
        output_xls_to_interval_negative_peaks_file = outputs['output_xls_to_interval_negative_peaks_file']


    cmdline = "%s --format='%s' --name='%s' --gsize='%s' --bw='%s' --mfold %s %s %s %s" % ( cmdline, options['format'], experiment_name, options['gsize'], options['bw'], options['mfoldlo'], options['mfoldhi'], options['nolambda'], options['bdg'] )


    if 'pvalue' in options:
        cmdline = "%s --pvalue='%s'" % ( cmdline, options['pvalue'] )
    elif 'qvalue' in options:
        cmdline = "%s --qvalue='%s'" % ( cmdline, options['qvalue'] )

    if 'broad' in options:
        cmdline = "%s --broad --broad-cutoff='%s'" % ( cmdline, options['broad_cutoff'] )

    if 'nomodel' in options:
        cmdline = "%s --nomodel --shiftsize='%s'" % ( cmdline, options['nomodel'] )
    #=================================================================================
    if options['command'] == "bdgcmp":
        output_bdgcmp = outputs['output_bdgcmp_file']

        cmdline = "%s -m %s -p %s -o bdgcmp_out.bdg" % ( cmdline, options['m'], options['pseudocount'] )
    #=================================================================================

    tmp_dir = tempfile.mkdtemp() #macs makes very messy output, need to contain it into a temp dir, then provide to user
    stderr_name = tempfile.NamedTemporaryFile().name # redirect stderr here, macs provides lots of info via stderr, make it into a report
    proc = subprocess.Popen( args=cmdline, shell=True, cwd=tmp_dir, stderr=open( stderr_name, 'wb' ) )
    proc.wait()
    #We don't want to set tool run to error state if only warnings or info, e.g. mfold could be decreased to improve model, but let user view macs log
    #Do not terminate if error code, allow dataset (e.g. log) creation and cleanup
    if proc.returncode:
        stderr_f = open( stderr_name )
        while True:
            chunk = stderr_f.read( CHUNK_SIZE )
            if not chunk:
                stderr_f.close()
                break
            sys.stderr.write( chunk )

    #=================================================================================
    #copy files created by macs2 to appripriate directory with the provided names
    #=================================================================================

    #=================================================================================
    #move files generated by callpeak command 
    if options['command'] == "callpeak":
        #run R to create pdf from model script
        if os.path.exists( os.path.join( tmp_dir, "%s_model.r" % experiment_name ) ):
            cmdline = 'R --vanilla --slave < "%s_model.r" > "%s_model.r.log"' % ( experiment_name, experiment_name )
            proc = subprocess.Popen( args=cmdline, shell=True, cwd=tmp_dir )
            proc.wait()

        #move bed out to proper output file
        created_bed_name =  os.path.join( tmp_dir, "%s_peaks.bed" % experiment_name )
        if os.path.exists( created_bed_name ):
            shutil.move( created_bed_name, output_bed )

        #OICR peak_xls file
        created_peak_xls_file =  os.path.join( tmp_dir, "%s_peaks.xls" % experiment_name )
        if os.path.exists( created_peak_xls_file ):
            shutil.copyfile( created_peak_xls_file, output_peaks )

        #peaks.encodepeaks (narrowpeaks) file
        created_narrowpeak_file = os.path.join (tmp_dir, "%s_peaks.encodePeak" % experiment_name )
        if os.path.exists( created_narrowpeak_file ):
            shutil.move (created_narrowpeak_file, output_narrowpeaks )

        #parse xls files to interval files as needed
        #if 'xls_to_interval' in options:
        if options['xls_to_interval'] == "True":
            create_peak_xls_file = os.path.join( tmp_dir, '%s_peaks.xls' % experiment_name )
            if os.path.exists( create_peak_xls_file ):
                xls_to_interval( create_peak_xls_file, output_xls_to_interval_peaks_file, header = 'peaks file' )
            create_peak_xls_file = os.path.join( tmp_dir, '%s_negative_peaks.xls' % experiment_name )
            if os.path.exists( create_peak_xls_file ):
                #print "negative file exists"
                xls_to_interval( create_peak_xls_file, output_xls_to_interval_negative_peaks_file, header = 'negative peaks file' )

        #move all remaining files to extra files path of html file output to allow user download
        out_html = open( output_extra_html, 'wb' )
        out_html.write( '<html><head><title>Additional output created by MACS (%s)</title></head><body><h3>Additional Files:</h3><p><ul>\n' % experiment_name )
        os.mkdir( output_extra_path )
        for filename in sorted( os.listdir( tmp_dir ) ):
            shutil.move( os.path.join( tmp_dir, filename ), os.path.join( output_extra_path, filename ) )
            out_html.write( '<li><a href="%s">%s</a></li>\n' % ( filename, filename ) )
        #out_html.write( '<li><a href="%s">%s</a>peakxls %s SomethingDifferent tmp_dir %s path %s exp_name %s</li>\n' % ( created_peak_xls_file, filename, filename, tmp_dir, output_extra_path, experiment_name ) )
        out_html.write( '</ul></p>\n' )
        out_html.write( '<h3>Messages from MACS:</h3>\n<p><pre>%s</pre></p>\n' % open( stderr_name, 'rb' ).read() )
        out_html.write( '</body></html>\n' )
        out_html.close()

    #=================================================================================
    #move files generated by bdgcmp command
    if options['command'] == "bdgcmp":
        created_bdgcmp_file = os.path.join (tmp_dir, "bdgcmp_out.bdg" )
        if os.path.exists( created_bdgcmp_file ):
            shutil.move (created_bdgcmp_file, output_bdgcmp )

    #=================================================================================    
    #cleanup
    #=================================================================================    
    os.unlink( stderr_name )
    os.rmdir( tmp_dir )

if __name__ == "__main__": main()
