#!/usr/bin/env python
import json
import optparse
import os
import re
import subprocess
import sys


def fetch_databases(genome_list=None):
    snpDBs = dict()
    databases_path = 'databases.out'
    args = ['snpEff', 'databases']
    with open(databases_path, 'w') as databases_output:
        return_code = subprocess.call(args=args, shell=False, stdout=databases_output.fileno())
    if return_code:
        sys.exit(return_code)
    try:
        with open(databases_path, 'r') as fh:
            for line in fh:
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
        sys.exit('Error parsing %s %s\n' % (databases_path, str(e)))
    return snpDBs


def getOrganismNames(genomes, organisms):
    genome_list = genomes.split(',')
    organism_list = organisms.split(',') if organisms else []
    if len(genome_list) != len(organism_list):
        descriptions = []
        snpDBdict = fetch_databases(genome_list=genome_list)
        for genome in snpDBdict:
            descriptions.append(snpDBdict[genome] if genome in snpDBdict else genome)
        return ','.join(descriptions)
    return organisms


def getSnpeffVersion():
    snpeff_version = 'SnpEff ?.?'
    stderr_path = 'snpeff.err'
    args = ['snpEff', '-h']
    with open(stderr_path, 'w') as stderr_fh:
        return_code = subprocess.call(args=args, shell=False, stderr=stderr_fh.fileno())
    if return_code != 255:
        sys.exit(return_code)
    with open(stderr_path) as fh:
        for line in fh:
            m = re.match(r'^[Ss]npEff version (SnpEff)\s*(\d+\.\d+).*$', line)
            if m:
                snpeff_version = m.groups()[0] + m.groups()[1]
                break
    return snpeff_version


# Download human database 'hg19'
# java -jar snpEff.jar download -v hg19
#
#        <command>java -jar \$SNPEFF_JAR_PATH/snpEff.jar download -c \$JAVA_JAR_PATH/snpEff.config $genomeVersion > $logfile </command>
#
# snpEffectPredictor.bin
# regulation_HeLa-S3.bin
# regulation_pattern = 'regulation_(.+).bin'
def download_database(data_manager_dict, target_directory, genome_version, organism):
    # get data_dir from config
    # ---
    # Databases are stored here
    # E.g.: Information for 'hg19' is stored in data_dir/hg19/
    #
    # Note: Since version 2.1 you can use tilde ('~') as first character to refer to your home directory
    data_dir = target_directory
    args = ['snpEff', 'download', '-dataDir', data_dir, '-v', genome_version]
    return_code = subprocess.call(args=args, shell=False)
    if return_code:
        sys.exit(return_code)
    # search data_dir/genome_version for files
    regulation_pattern = 'regulation_(.+).bin'
    genome_path = os.path.join(data_dir, genome_version)
    snpeff_version = getSnpeffVersion()
    key = snpeff_version + '_' + genome_version
    if os.path.isdir(genome_path):
        for _, _, files in os.walk(genome_path):
            for fname in files:
                if fname.startswith('snpEffectPredictor'):
                    # if snpEffectPredictor.bin download succeeded
                    name = genome_version + (' : ' + organism if organism else '')
                    data_table_entry = dict(key=key, version=snpeff_version, value=genome_version, name=name, path=data_dir)
                    _add_data_table_entry(data_manager_dict, 'snpeffv_genomedb', data_table_entry)
                else:
                    m = re.match(regulation_pattern, fname)
                    if m:
                        name = m.groups()[0]
                        data_table_entry = dict(key=key, version=snpeff_version, genome=genome_version, value=name, name=name)
                        _add_data_table_entry(data_manager_dict, 'snpeffv_regulationdb', data_table_entry)
    return data_manager_dict


def _add_data_table_entry(data_manager_dict, data_table, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table] = data_manager_dict['data_tables'].get(data_table, [])
    data_manager_dict['data_tables'][data_table].append(data_table_entry)
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    parser.add_option('-g', '--genome_version', dest='genome_version', action='store', type="string", default=None, help='genome_version')
    parser.add_option('-o', '--organism', dest='organism', action='store', type="string", default=None, help='organism name')
    (options, args) = parser.parse_args()

    filename = args[0]

    with open(filename) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    data_manager_dict = {}

    # Create SnpEff Reference Data
    for genome_version, organism in zip(options.genome_version.split(','), getOrganismNames(options.genome_version, options.organism).split(',')):
        download_database(data_manager_dict, target_directory, genome_version, organism)

    # save info to json file
    with open(filename, 'w'):
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
