#!/usr/bin/env python

import argparse
import glob
import json
import logging
import os


import eutils


logging.basicConfig(level=logging.INFO)


def handleEfetchException(e, db, payload):
    logging.error('No results returned. This could either be due to no records matching the supplied IDs for the query database or it could be an error due to invalid parameters.  The reported exception was "%s".\n\nPayload used for the efetch query to database "%s"\n\n%s', e, db, json.dumps(payload, indent=4))

    # Create a file in the downloads folder so that the user can access run information
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'downloads')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)

    print('The following files were downloaded:')
    print(os.listdir(final_directory))

    file_path = os.path.join('downloads', 'no_results.txt')
    with open(file_path, 'w') as handle:
        handle.write('No results')


def localFetch(db, gformat, newname, **payload):
    problem = None
    try:
        c.fetch(db, **payload)

        for chunk, file in enumerate(glob.glob('downloads/EFetch *')):
            os.rename(file, '%s%s.%s' % (newname, chunk + 1, gformat))

    except Exception as e:
        problem = e
        handleEfetchException(e, db, payload)
    else:
        print('The following files were downloaded:')
        print(os.listdir('downloads'))

    return problem


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EFetch', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    parser.add_argument('--version', action='version', version=eutils.Client.getVersion(), help='Version (reports Biopython version)')

    # ID source
    parser.add_argument('--id_json', help='list of ids in a json file as returned by esearch or elink')
    parser.add_argument('--id_xml', help='list of ids in an xml file as returned by esearch or elink')
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Fetch results from previous query (JSON)')
    parser.add_argument('--history_xml', help='Fetch results from previous query (XML)')

    # Output
    parser.add_argument('--retmode', help='Retmode')
    parser.add_argument('--rettype', help='Rettype')
    parser.add_argument('--galaxy_format', help='Galaxy format')
    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)

    payload = {}
    for attr in ('retmode', 'rettype'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    if args.history_file is not None or args.history_xml is not None:
        if args.history_file is not None:
            input_histories = c.get_histories()
        else:
            input_histories = c.extract_histories_from_xml_file(args.history_xml)

        problem = None
        for hist in input_histories:
            qkey = hist['query_key']
            tmp_payload = payload
            tmp_payload.update(hist)
            newname = 'downloads/EFetch-%s-%s-querykey%s-chunk' % (args.rettype, args.retmode, qkey)
            problem = localFetch(args.db, args.galaxy_format, newname, **tmp_payload)

            if os.path.exists('downloads'):
                os.rename('downloads', 'downloads-qkey%s' % (qkey))

        if not os.path.exists('downloads'):
            os.makedirs('downloads')

        for relpath in glob.glob('downloads-qkey*/*'):
            file = os.path.basename(relpath)
            os.rename(relpath, 'downloads/%s' % (file))

        if problem is not None:
            raise problem

    else:
        merged_ids = c.parse_ids(args.id_list, args.id, args.history_file, args.id_xml, args.id_json)
        payload['id'] = ','.join(merged_ids)
        newname = 'downloads/EFetch-%s-%s-chunk' % (args.rettype, args.retmode)
        localFetch(args.db, args.galaxy_format, newname, **payload)
