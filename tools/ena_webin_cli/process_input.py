import gzip
import json
import os
import sys
import shutil
import yaml

def get_section_string(f, start_line, end_line, return_string=False):
    # consume starting lines
    start_string = iter(f.readline, start_line)
    start_string = ''.join(line for line in start_string)
    # read YAML lines
    yaml_string = iter(f.readline, end_line)
    if return_string:
        return ''.join(x for x in yaml_string)
    else:
        return [x for x in yaml_string]
  
def fill_from_yaml_data(yaml_only_dict, studies_samples_dict):
    # fill experiment information (platform)  **** 
    for index,exp in yaml_only_dict['ENA_experiment'].items():
        study_alias = exp['study_alias']
        sample_alias = exp['sample_alias']
        if study_alias in studies_samples_dict.keys():
            if sample_alias in studies_samples_dict[study_alias].keys():
                studies_samples_dict[study_alias][sample_alias]['experiments'].append({'platform': exp['platform']})
            else:
                studies_samples_dict[study_alias][sample_alias] = {'experiments': [{'platform': exp['platform']}]}
        else:
            studies_samples_dict[study_alias] = {sample_alias: {'experiments':[{'platform': exp['platform']}]}}


def load_receipt_data(input_file_path):
    # should do some health check of the input file?
    # load yaml section
    loaded_data = {} 
    yaml_delimiter = 'YAML -------------\n'
    with open(input_file_path) as input_file:
        yaml_only_section = yaml.safe_load(get_section_string(input_file, start_line=yaml_delimiter, end_line=yaml_delimiter, return_string=True))
    fill_from_yaml_data(yaml_only_section, loaded_data)
    # read study accessions
    study_delimiter = 'Study accession details:\n'
    end_line = '\n'
    with open(input_file_path) as input_file:
        studies_accession_lines = get_section_string(input_file, start_line=study_delimiter, end_line=end_line)
    # loaded_data['studies'] = {}
    for study_line in studies_accession_lines:
        if study_line != '\n':
            alias, accession, *_ = study_line.split('\t')
            try:
                loaded_data[alias]['accession'] = accession
            except KeyError:
                print(f"Experiment {exp} has unknown study or sample")
            # loaded_data['studies'][alias]['accession'] = accession
    samples_delimiter = 'Sample accession details:\n'
    with open(input_file_path) as input_file:
        samples_accession_lines = get_section_string(input_file, start_line=samples_delimiter, end_line=end_line)
        ## need to iterate over all studies, because here I don't know which study is the sample from.
    # loaded_data['samples'] = {}
    for sample_line in samples_accession_lines:
        if sample_line != '\n':
            alias, accession, *_ = sample_line.split('\t')
            for study in loaded_data.keys():
                if alias in loaded_data[study].keys():
                    loaded_data[study][alias]['accession'] = accession
                    break
    return loaded_data


"""
Takes as input:
    1. A receipt obtained from ENA submission tool: 
        a txt file that contains sections describing submission details.
    2. A json file with the list of fasta that the user loaded
    3. Path to write generated manifests
    4. Manifest template path: the manifest with the global values set 
        (e.g COVERAGE, MINGAPLENGHT..)
"""



def main():
    input_file_path = sys.argv[1]
    fasta_names_list_path = sys.argv[2]
    out_manifest_base = sys.argv[3] 
    manifest_template = sys.argv[4]
    # load submitted data from receipt file
    data_dict = load_receipt_data(input_file_path)
    # iterate over the list of fasta files
    with open(fasta_names_list_path, 'r') as fasta_files_json_file:
        fasta_files_list = json.load(fasta_files_json_file)
    with open('submit_list.tab', 'w') as written_manifests_out:
        for fasta_file in fasta_files_list:
            if fasta_file.endswith('.fasta.gz'):
                sample_alias = fasta_file[:-9]
            else:
                sample_alias = fasta_file[:-6]
            print(f'Processing {sample_alias}')
            found_metadata = False
            for study_alias in data_dict.keys():
                if sample_alias in data_dict[study_alias].keys():
                    sample_accession = data_dict[study_alias][sample_alias]['accession']
                    study_accession = data_dict[study_alias]['accession']
                    ### TODO get a string that concatenates plaform information from multiple exp
                    platform = data_dict[study_alias][sample_alias]['experiments'][0]['platform']
                    manifest_path = os.path.join(out_manifest_base, sample_alias + '.manifest.txt')
                    with open(manifest_path, "w") as output_handle:
                        # first dump the contents of manifest template
                        # containing the global vars
                        with open(manifest_template) as m_template:
                            output_handle.write(m_template.read())
                        output_handle.write("ASSEMBLYNAME\tconsensus_" + sample_alias + "\n")
                        output_handle.write("PLATFORM\t" + platform + "\n")
                        output_handle.write("STUDY\t" + study_accession + "\n")
                        output_handle.write("SAMPLE\t" + sample_accession + "\n")
                        # files should be available in the corresponding dir and named:
                        #  sample_alias.fasta.gz 
                        output_handle.write("FASTA\t" + sample_alias + '.fasta.gz' + "\n")
                    found_metadata = True
                    written_manifests_out.write(manifest_path + '\n')
                    break
            if not found_metadata:
                print(f'No metadata found for sample {sample_alias}')


if __name__ == '__main__':
    main()
