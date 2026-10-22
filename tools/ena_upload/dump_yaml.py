import sys

import yaml


def fetch_table_data(table_path):
    data_dict = {}
    with open(table_path) as table_to_load:
        # load headers
        headers = table_to_load.readline().strip('\n').split('\t')
        row_id = 0
        for line in table_to_load.readlines():
            # print(line)
            line_data = line.strip('\n').split('\t')
            row_dict = {}
            for col_num in range(len(headers)):
                col_name = headers[col_num]
                row_dict[col_name] = line_data[col_num]
            data_dict[row_id] = row_dict
            row_id += 1
        return data_dict


all_data_dict = {}
print('YAML -------------')
studies_table_path = sys.argv[1]
table_data = fetch_table_data(studies_table_path)
all_data_dict['ENA_study'] = table_data
samples_table_path = sys.argv[2]
table_data = fetch_table_data(samples_table_path)
all_data_dict['ENA_sample'] = table_data
experiments_table_path = sys.argv[3]
table_data = fetch_table_data(experiments_table_path)
all_data_dict['ENA_experiment'] = table_data
runs_table_path = sys.argv[4]
table_data = fetch_table_data(runs_table_path)
all_data_dict['ENA_run'] = table_data
# print(all_data_dict)
print(yaml.dump(all_data_dict))
print('YAML -------------')
