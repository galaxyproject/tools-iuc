#!/usr/bin/env python
# Dan Blankenberg


import argparse
import json
import os
import re
import shutil
import urllib.request


OUTPUT_BASE_NAME = 'RAVmodel.rds'

BUCKET_NAME = "genomic_super_signature"
RAVMODEL_NAME = "RAVmodel_{prior}.rds"
PRIOR_URL = "https://storage.googleapis.com/{bucket_name}/{prior_filename}"


def get_prior_url(prior):
    return PRIOR_URL.format(bucket_name=BUCKET_NAME, prior_filename=RAVMODEL_NAME.format(prior=prior))


def clean(text):
    return re.sub(r'[^\w\-_\.]', '-', text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', dest='output', action='store', type=str, default=None)
    parser.add_argument('--json', dest='json', action='store', type=str, default=None)
    parser.add_argument('--dbkey', dest='dbkey', action='store', type=str, default=None, help='dbkey')
    parser.add_argument('--name', dest='name', action='store', type=str, default=None)
    parser.add_argument('--id', dest='id', action='store', type=str, default=None)
    parser.add_argument('--prior', dest='prior', action='store', type=str, default=None)
    parser.add_argument('--url', dest='url', action='store', type=str, default=None)
    parser.add_argument('--input', dest='input', action='store', type=str, default=None)
    parser.add_argument('--symlink', dest='symlink', action='store_true', default=False)
    parser.add_argument('--version', dest='version', action='store', type=str, default="0")
    parser.add_argument('--tablename', dest='tablename', action='store', type=str, default="genomic_super_signature_ravmodels")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    output_filename = os.path.join(args.output, OUTPUT_BASE_NAME)
    if args.symlink:
        assert args.input, "You must provide a filename when using symlink functionality."
        os.symlink(args.input, output_filename)
    else:
        url = args.url
        if args.prior:
            url = get_prior_url(args.prior)
        if url:
            urllib.request.urlretrieve(url, output_filename)
        else:
            assert args.input, "You must provide a filename, prior, or URL."
            shutil.copyfile(args.input, output_filename)
    data_manager_dict = {'data_tables': {args.tablename: [dict(value=clean(args.id), dbkey=args.dbkey, version=args.version, name=args.name, path=OUTPUT_BASE_NAME)]}}

    with open(args.json, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
