#!/usr/bin/env python

import sys
import os
import re
import tempfile
import subprocess
import fileinput
import shutil
import optparse
import urllib2
from ftplib import FTP
import tarfile

from galaxy.util.json import from_json_string, to_json_string

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)

"""
# Download human database 'hg19'
java -jar snpEff.jar download -v hg19

        <command>java -jar \$SNPEFF_JAR_PATH/snpEff.jar download -c \$JAVA_JAR_PATH/snpEff.config $genomeVersion > $logfile </command>

snpEffectPredictor.bin
regulation_HeLa-S3.bin
regulation_pattern = 'regulation_(.+).bin'


"""
def download_database(data_manager_dict, target_directory, jar_path,config,genome_version,organism):
    ## get data_dir from config 
    ##---
    ## Databases are stored here
    ## E.g.: Information for 'hg19' is stored in data_dir/hg19/
    ##
    ## Note: Since version 2.1 you can use tilde ('~') as first character to refer to your home directory
    ##---
    #data_dir = ~/snpEff/data/
    data_dir = target_directory
    (snpEff_dir,snpEff_jar) = os.path.split(jar_path)
    args = [ 'java','-jar' ]
    args.append( jar_path )
    args.append( 'download' )
    args.append( '-c' )
    args.append( config )
    args.append( '-dataDir' )
    args.append( data_dir )
    args.append( '-v' )
    args.append( genome_version )
    proc = subprocess.Popen( args=args, shell=False, cwd=snpEff_dir )
    return_code = proc.wait()
    if return_code:
        sys.exit( return_code )
    ## search data_dir/genome_version for files
    regulation_pattern = 'regulation_(.+).bin'
    #  annotation files that are included in snpEff by a flag
    annotations_dict = {'nextProt.bin' : '-nextprot','motif.bin': '-motif'}
    genome_path = os.path.join(data_dir,genome_version)
    if os.path.isdir(genome_path):
        for root, dirs, files in os.walk(genome_path):
            for fname in files:
                if fname.startswith('snpEffectPredictor'):
                    # if snpEffectPredictor.bin download succeeded
                    name = genome_version + (' : ' + organism if organism else '') 
                    data_table_entry = dict(value=genome_version, name=name, path=data_dir)
                    _add_data_table_entry( data_manager_dict, 'snpeff_genomedb', data_table_entry )
                else:
                    m = re.match(regulation_pattern,fname)
                    if m:
                        name = m.groups()[0]
                        data_table_entry = dict(genome=genome_version,value=name, name=name)
                        _add_data_table_entry( data_manager_dict, 'snpeff_regulationdb', data_table_entry )
                    elif fname in annotations_dict:
                        value = annotations_dict[fname]
                        name = value.lstrip('-')
                        data_table_entry = dict(genome=genome_version,value=value, name=name)
                        _add_data_table_entry( data_manager_dict, 'snpeff_annotations', data_table_entry )
    return data_manager_dict

def _add_data_table_entry( data_manager_dict, data_table, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table] = data_manager_dict['data_tables'].get( data_table, [] )
    data_manager_dict['data_tables'][data_table].append( data_table_entry )
    return data_manager_dict

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-j', '--jar_path', dest='jar_path', action='store', type="string", default=None, help='snpEff.jar path' )
    parser.add_option( '-c', '--config', dest='config', action='store', type="string", default=None, help='snpEff.config path' )
    parser.add_option( '-g', '--genome_version', dest='genome_version', action='store', type="string", default=None, help='genome_version' )
    parser.add_option( '-o', '--organism', dest='organism', action='store', type="string", default=None, help='organism name' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = from_json_string( open( filename ).read() )
    target_directory = params[ 'output_data' ][0][.files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}


    #Create SnpEff Reference Data
    for genome_version, organism in zip(options.genome_version.split(','), options.organism.split(',')):
        download_database( data_manager_dict, target_directory, options.jar_path, options.config, genome_version, organism )

    #save info to json file
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )

if __name__ == "__main__": main()

