#!/usr/bin/env python

#from __future__ import print_function

import json
import optparse

def _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get(data_table_name, [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    parser.add_option('-o', '--organism', dest='organism', action='store', type="string", default=None, help='organism')
    parser.add_option('-n', '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name')
    parser.add_option('-d', '--DB_path', dest='DB_path', action='store', type="string", default=None, help='DB_path')
    parser.add_option('--organism_version', dest='organism_version', action='store', type="string", default=None, help='organism_version')
    parser.add_option('--promoter_version', dest='promoter_version', action='store', type="string", default=None, help='promoter_version')
    (options, args) = parser.parse_args()

    filename = args[0]
    with open(filename) as fh:
        params = json.load(fh)
    
    dbkey = str(options.organism) + '_o' + str(options.organism_version) + '_p' + str(options.promoter_version)
    data_manager_dict = {}
    data_table_entry = dict(value=dbkey, dbkey=dbkey, organism=options.organism, path=options.DB_path, organism_version=options.organism_version, promoter_version=options.promoter_version)
    _add_data_table_entry(data_manager_dict, options.data_table_name, data_table_entry)

    # Save info to json file
    with open(filename, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
