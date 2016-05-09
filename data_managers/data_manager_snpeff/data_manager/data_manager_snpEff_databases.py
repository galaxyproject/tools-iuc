#!/usr/bin/env python
import json
import optparse
import os
import subprocess
import sys


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def fetch_databases(data_manager_dict, target_directory, jar_path):
    (snpEff_dir, snpEff_jar) = os.path.split(jar_path)
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    databases_path = os.path.join( target_directory, 'databases.out' )
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
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['snpeff4_databases'] = data_manager_dict['data_tables'].get( 'snpeff4_databases', [] )
    data_table_entries = []
    try:
        fh = open(databases_path, 'r')
        for i, line in enumerate(fh):
            fields = line.split('\t')
            if len(fields) >= 2:
                genome_version = fields[0].strip()
                if genome_version.startswith("Genome") or genome_version.startswith("-"):
                    continue
                # snpeff test genome
                if genome_version == '30c2c903' or fields[1].strip() == 'TestCase' or fields[1].strip().startswith('Test_'):
                    continue
                description = fields[1].strip() + ' : ' + genome_version
                data_table_entries.append(dict(value=genome_version, name=description))
        data_manager_dict['data_tables']['snpeff4_databases'] = data_table_entries
    except Exception as e:
        stop_err( 'Error parsing %s %s\n' % (databases_path, str( e )) )
    else:
        fh.close()
    return data_manager_dict


def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-j', '--jar_path', dest='jar_path', action='store', type="string", default=None, help='snpEff.jar path' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = json.loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    # Create Defuse Reference Data
    data_manager_dict = fetch_databases( data_manager_dict, target_directory, options.jar_path)

    # save info to json file
    open( filename, 'wb' ).write( json.dumps( data_manager_dict ) )

if __name__ == "__main__":
    main()
