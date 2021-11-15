import json

import requests

URL = "https://www.ebi.ac.uk/ena/portal/api/search"


def check_remote_entry(entry_type, query_dict, out_format='json'):
    '''
    Checks if an entry with that alias exists in the ENA repos
    entry_type = [study | sample | experiment | run]
    '''
    assert entry_type in ['study', 'sample', 'experiment', 'run']
    params_dict = {}
    query_str = ' AND '.join(['%s="%s"' % (key, value) for (key, value) in query_dict.items()])
    params_dict['query'] = query_str
    params_dict['result'] = 'read_' + entry_type
    params_dict['fields'] = entry_type + '_alias'
    params_dict['format'] = out_format
    response = requests.post(URL, data=params_dict)
    if response.content != b'':
        return json.loads(response.content)
    return []
