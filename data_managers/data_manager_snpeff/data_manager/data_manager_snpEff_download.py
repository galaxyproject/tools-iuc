#!/usr/bin/env python
import gzip
import json
import optparse
import os
import re
import subprocess
import sys


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def fetch_databases(jar_path, genome_list=None):
    snpDBs = dict()
    (snpEff_dir, snpEff_jar) = os.path.split(jar_path)
    databases_path = 'databases.out'
    databases_output = open(databases_path, 'w')
    args = [ 'java', '-jar' ]
    args.append( snpEff_jar )
    args.append( 'databases' )
    # tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-data-manager-snpEff-stderr" )
    # databases_output = open(databases_path)
    # proc = subprocess.Popen( args=args, shell=False, cwd=snpEff_dir, stdout=databases_output.fileno(), stderr=tmp_stderr.fileno() )
    proc = subprocess.Popen( args=args, shell=False, cwd=snpEff_dir, stdout=databases_output.fileno() )
    return_code = proc.wait()
    if return_code:
        sys.exit( return_code )
    databases_output.close()
    try:
        fh = open(databases_path, 'r')
        for i, line in enumerate(fh):
            fields = line.split('\t')
            if len(fields) >= 2:
                genome_version = fields[0].strip()
                if genome_list and genome_version not in genome_list:
                    continue
                if genome_version.startswith("Genome") or genome_version.startswith("-"):
                    continue
                description = fields[1].strip()
                snpDBs[genome_version] = description
    except Exception as e:
        stop_err( 'Error parsing %s %s\n' % (databases_path, str( e )) )
    else:
        fh.close()
    return snpDBs


def getOrganismNames(jar_path, genomes, organisms):
    genome_list = genomes.split(',')
    organism_list = organisms.split(',') if organisms else []
    if len(genome_list) != len(organism_list):
        descriptions = []
        snpDBdict = fetch_databases(jar_path, genome_list=genome_list)
        for genome in snpDBdict:
            descriptions.append(snpDBdict[genome] if genome in snpDBdict else genome)
        return ','.join(descriptions)
    return organisms


def getSnpeffVersion(jar_path):
    snpeff_version = 'SnpEff ?.?'
    (snpEff_dir, snpEff_jar) = os.path.split(jar_path)
    stderr_path = 'snpeff.err'
    stderr_fh = open(stderr_path, 'w')
    args = [ 'java', '-jar' ]
    args.append( snpEff_jar )
    args.append( '-h' )
    proc = subprocess.Popen( args=args, shell=False, cwd=snpEff_dir, stderr=stderr_fh.fileno() )
    return_code = proc.wait()
    if return_code != 255:
        sys.exit( return_code )
    stderr_fh.close()
    fh = open(stderr_path, 'r')
    for line in fh:
        m = re.match('^[Ss]npEff version (SnpEff)\s*(\d+\.\d+).*$', line)
        if m:
            snpeff_version = m.groups()[0] + m.groups()[1]
            break
    fh.close()
    return snpeff_version


# Starting with SnpEff 4.1 the .bin files contain the SnpEff version:
# Example - the first 3 line of GRCh37.75/snpEffectPredictor.bin (uncompressed):
#
# SnpEff  4.1
# CHROMOSOME      2       1       0       179197  GL000219.1      false
# CHROMOSOME      3       1       0       81347269        HSCHR17_1       false
def getSnpeffVersionFromFile(path):
    snpeff_version = None
    try:
        fh = gzip.open(path, 'rb')
        buf = fh.read(100)
        lines = buf.splitlines()
        m = re.match('^(SnpEff)\s+(\d+\.\d+).*$', lines[0].strip())
        if m:
            snpeff_version = m.groups()[0] + m.groups()[1]
        fh.close()
    except Exception as e:
        stop_err( 'Error parsing SnpEff version from: %s %s\n' % (path, str( e )) )
    return snpeff_version


# Download human database 'hg19'
# java -jar snpEff.jar download -v hg19
#
#        <command>java -jar \$SNPEFF_JAR_PATH/snpEff.jar download -c \$JAVA_JAR_PATH/snpEff.config $genomeVersion > $logfile </command>
#
# snpEffectPredictor.bin
# regulation_HeLa-S3.bin
# regulation_pattern = 'regulation_(.+).bin'
def download_database(data_manager_dict, target_directory, jar_path, config, genome_version, organism):
    # get data_dir from config
    # ---
    # Databases are stored here
    # E.g.: Information for 'hg19' is stored in data_dir/hg19/
    #
    # Note: Since version 2.1 you can use tilde ('~') as first character to refer to your home directory
    # ---
    # data_dir = ~/snpEff/data/
    data_dir = target_directory
    (snpEff_dir, snpEff_jar) = os.path.split(jar_path)
    args = [ 'java', '-jar' ]
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
    # search data_dir/genome_version for files
    regulation_pattern = 'regulation_(.+).bin'
    #  annotation files that are included in snpEff by a flag
    annotations_dict = {'nextProt.bin': '-nextprot', 'motif.bin': '-motif'}
    genome_path = os.path.join(data_dir, genome_version)
    snpeff_version = getSnpeffVersion(jar_path)
    key = snpeff_version + '_' + genome_version
    if os.path.isdir(genome_path):
        for root, dirs, files in os.walk(genome_path):
            for fname in files:
                if fname.startswith('snpEffectPredictor'):
                    # if snpEffectPredictor.bin download succeeded
                    name = genome_version + (' : ' + organism if organism else '')
                    # version = getSnpeffVersionFromFile(os.path.join(root,fname))
                    data_table_entry = dict(key=key, version=snpeff_version, value=genome_version, name=name, path=data_dir)
                    _add_data_table_entry( data_manager_dict, 'snpeffv_genomedb', data_table_entry )
                else:
                    m = re.match(regulation_pattern, fname)
                    if m:
                        name = m.groups()[0]
                        data_table_entry = dict(key=key, version=snpeff_version, genome=genome_version, value=name, name=name)
                        _add_data_table_entry( data_manager_dict, 'snpeffv_regulationdb', data_table_entry )
                    elif fname in annotations_dict:
                        value = annotations_dict[fname]
                        name = value.lstrip('-')
                        data_table_entry = dict(key=key, version=snpeff_version, genome=genome_version, value=value, name=name)
                        _add_data_table_entry( data_manager_dict, 'snpeffv_annotations', data_table_entry )
    return data_manager_dict


def _add_data_table_entry( data_manager_dict, data_table, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table] = data_manager_dict['data_tables'].get( data_table, [] )
    data_manager_dict['data_tables'][data_table].append( data_table_entry )
    return data_manager_dict


def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-j', '--jar_path', dest='jar_path', action='store', type="string", default=None, help='snpEff.jar path' )
    parser.add_option( '-c', '--config', dest='config', action='store', type="string", default=None, help='snpEff.config path' )
    parser.add_option( '-g', '--genome_version', dest='genome_version', action='store', type="string", default=None, help='genome_version' )
    parser.add_option( '-o', '--organism', dest='organism', action='store', type="string", default=None, help='organism name' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    # Create SnpEff Reference Data
    for genome_version, organism in zip(options.genome_version.split(','), getOrganismNames(options.jar_path, options.genome_version, options.organism).split(',')):
        download_database( data_manager_dict, target_directory, options.jar_path, options.config, genome_version, organism )

    # save info to json file
    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
