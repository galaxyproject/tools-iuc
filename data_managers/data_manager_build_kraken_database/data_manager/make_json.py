import argparse
import json
import os
import shutil


def main(args):
    data_manager_entry = {}
    data_manager_entry['value'] = args.database.lower()
    data_manager_entry['name'] = args.database
    data_manager_entry['path'] = '.'
    data_manager_json = dict(data_tables=dict(kraken_databases=data_manager_entry))
    params = json.loads(open(args.output).read())
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    output_path = os.path.join(os.getcwd(), 'kraken-database')
    for filename in os.listdir(output_path):
        shutil.move(os.path.join(output_path, filename), target_directory)
    open(args.output, 'w').write(json.dumps(data_manager_json))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--db', dest='database', action='store', help='Database name')
    parser.add_argument('--out', dest='output', action='store', help='JSON filename')
    args = parser.parse_args()
    main(args)
