#!/usr/bin/env python
import gzip
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
    args = ['snpEff', '-version']
    try:
        version_output = subprocess.check_output(args, shell=False).decode()
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    m = re.match(r'^(SnpEff)\s*(\d+\.\d+).*$', version_output)
    if m:
        snpeff_version = m.groups()[0] + m.groups()[1]
    return snpeff_version


def getSnpeffDbVersion(db_predictor_bin_file):
    snpeff_db_version = None
    try:
        with gzip.open(db_predictor_bin_file, "rt") as fh:
            buf = fh.read(100)
            lines = buf.splitlines()
            m = re.match(r"^(SnpEff)\s+(\d+\.\d+).*$", lines[0].strip())
            if m:
                snpeff_db_version = m.groups()[0] + m.groups()[1]
    except Exception:
        pass
    return snpeff_db_version


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
    db_version = None
    genomedb_name = regulationdb_name = ""

    if os.path.isdir(genome_path):
        for dirpath, _, files in os.walk(genome_path):
            for fname in files:
                if fname.startswith('snpEffectPredictor'):
                    # if snpEffectPredictor.bin download succeeded
                    genomedb_name = genome_version + (' : ' + organism if organism else '')
                    db_version = getSnpeffDbVersion(os.path.join(dirpath, fname)) or snpeff_version
                else:
                    m = re.match(regulation_pattern, fname)
                    if m:
                        regulationdb_name = m.groups()[0]

        if db_version:
            data_table_entry = dict(
                key=key,
                version=db_version,
                value=genome_version,
                name=genomedb_name,
                path=f"snpEff/{db_version}/data"
            )
            data_manager_dict['data_tables']['snpeffv_genomedb'].append(data_table_entry)

        if regulationdb_name:
            data_table_entry = dict(
                key=key,
                version=db_version or snpeff_version,
                genome=genome_version,
                value=regulationdb_name,
                name=regulationdb_name
            )
            data_manager_dict['data_tables']['snpeffv_regulationdb'].append(data_table_entry)


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
    data_manager_dict = {
        'data_tables': {
            'snpeffv_genomedb': [],
            'snpeffv_regulationdb': []
        }
    }

    # Create SnpEff Reference Data
    for genome_version, organism in zip(options.genome_version.split(','), getOrganismNames(options.genome_version, options.organism).split(',')):
        download_database(data_manager_dict, target_directory, genome_version, organism)

    # save info to json file
    with open(filename, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
