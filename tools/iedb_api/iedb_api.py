#!/usr/bin/env python

import argparse
import os.path
import re
import sys
import time
from urllib.error import HTTPError
from urllib.parse import unquote, urlencode
from urllib.request import urlopen

mhci_methods = ['recommended', 'consensus',
                'netmhcpan_ba', 'netmhcpan_el',
                'ann', 'smmpmbec', 'smm',
                'comblib_sidney2008', 'netmhccons',
                'pickpocket', 'netmhcstabpan']
mhcii_methods = ['recommended', 'consensus', 'NetMHCIIpan',
                 'nn_align', 'smm_align', 'comblib', 'tepitope']
processing_methods = ['recommended', 'netmhcpan', 'ann',
                      'smmpmbec', 'smm', 'comblib_sidney2008',
                      'netmhccons', 'pickpocket']
mhcnp_methods = ['mhcnp', 'netmhcpan']
bcell_methods = ['Bepipred', 'Chou-Fasman', 'Emini', 'Karplus-Schulz',
                 'Kolaskar-Tongaonkar', 'Parker', 'Bepipred-2.0']
prediction_methods = {'mhci': mhci_methods,
                      'mhcii': mhcii_methods,
                      'processing': processing_methods,
                      'mhcnp': mhcnp_methods,
                      'bcell': bcell_methods}
all_methods = set(mhci_methods + mhcii_methods + mhcnp_methods + bcell_methods)
prediction_lengths = {'mhci': range(8, 16),
                      'mhcii': range(11, 31),
                      'processing': range(8, 15),
                      'mhcnp': range(8, 12),
                      'bcell': range(8, 16)}


def parse_alleles(allelefile, query_lengths):
    alleles = []
    lengths = []
    with open(allelefile, 'r') as fh:
        for i, line in enumerate(fh):
            fields = line.strip().split(',')
            allele = fields[0].strip()
            if allele:
                if len(fields) > 1:
                    for alen in fields[1:]:
                        alleles.append(allele)
                        lengths.append(alen)
                elif query_lengths:
                    lens = []
                    for ql in query_lengths:
                        lens.extend(str(ql).split(','))
                    for alen in lens:
                        alleles.append(allele)
                        lengths.append(alen)
                else:
                    alleles.append(allele)
    return (alleles, lengths)


def query(url, prediction, seq, allele, length, results,
          seqid=None, method='recommended', proteasome=None,
          timeout=300, retries=3, sleep=300, debug=False):
    params = dict()
    if method:
        params['method'] = method.encode()
    if proteasome:
        params['proteasome'] = proteasome.encode()
    params['sequence_text'] = seq.strip().encode()
    if allele is not None:
        params['allele'] = allele.encode()
    if length is not None:
        if prediction == 'bcell':
            params['window_size'] = str(length).encode()
        elif length == 'asis':
            params['length'] = str(length).encode()
        else:
            slen = len(seq)
            alleles = []
            lengths = []
            for i in zip(length.split(','), allele.split(',')):
                if int(i[0]) <= slen:
                    lengths.append(i[0])
                    alleles.append(i[1])
            if lengths:
                params['length'] = str(','.join(lengths)).encode()
                params['allele'] = str(','.join(alleles)).encode()
            else:
                return results
    req_data = urlencode(params)
    if debug:
        print('url %s %s' % (url, unquote(req_data)), file=sys.stderr)
    retries = max(0, retries) + 1
    for retry in range(1, retries):
        response = None
        try:
            response = urlopen(url, data=req_data.encode('utf-8'),
                               timeout=timeout)
            if response and response.getcode() == 200:
                data = [line.decode() for line in response.readlines()]
                if debug:
                    print(data, file=sys.stderr)
                rslts = results['prediction']['entries']
                for ln, line in enumerate(data):
                    if 'invalid' in line.lower() or 'tools_api.html' in line:
                        msg = '%s %s\n%s' % (url, unquote(req_data),
                                             ''.join(data))
                        warn_err(msg, exit_code=1)
                    if line.find('eptide') > 0:
                        results['prediction']['header'] = "#%s%s" %\
                            ("ID\t" if seqid else "", line)
                        continue
                    elif method == 'Bepipred' and line.find('Residue') > 0:
                        results['detail']['header'] = "#%s%s" %\
                            ("ID\t" if seqid else "", line)
                        rslts = results['detail']['entries']
                        continue
                    if seqid:
                        rslts.extend("%s\t%s" % (seqid, line))
                    else:
                        rslts.extend(line)
                break
            else:
                code = response.getcode() if response else 1
                warn_err("Error connecting to IEDB server\n",
                         exit_code=code)
        except HTTPError as e:
            code = None if retry < retries else e.code
            warn_err("%d of %d Error connecting to IEDB server %s\n" %
                     (retry, retries, e),
                     exit_code=code)
            time.sleep(sleep)
        except Exception as e:
            warn_err("Error connecting to IEDB server %s\n" % e,
                     exit_code=3)
    return results


def warn_err(msg, exit_code=1):
    sys.stderr.write(msg)
    sys.stderr.flush()
    if exit_code:
        sys.exit(exit_code)


def __main__():
    # Parse Command Line
    parser = argparse.ArgumentParser(description='', epilog='')
    parser.add_argument('-p', '--prediction',
                        default='mhci',
                        choices=prediction_methods.keys(),
                        help='IEDB API prediction service')
    parser.add_argument('-s', '--sequence',
                        action="append",
                        default=None,
                        help='Peptide Sequence')
    parser.add_argument('-m', '--method',
                        default='recommended',
                        choices=all_methods,
                        help='prediction method')
    parser.add_argument('-P', '--proteasome',
                        default=None,
                        choices=['immuno', 'constitutive'],
                        help='IEDB processing proteasome type')
    parser.add_argument('-a', '--allele',
                        action="append",
                        default=[],
                        help='Alleles for which to make predictions')
    parser.add_argument('-A', '--allelefile',
                        default=None,
                        help='File of HLA alleles')
    parser.add_argument('-l', '--length',
                        action="append",
                        default=[],
                        help='lengths for which to make predictions for alleles')
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=None,
                        help='window_size for bcell prediction')
    parser.add_argument('-i', '--input',
                        default=None,
                        help='Input file for peptide sequences '
                             + '(fasta or tabular)')
    parser.add_argument('-c', '--column',
                        default=None,
                        help='Peptide Column in a tabular input file')
    parser.add_argument('-C', '--id_column',
                        default=None,
                        help='ID Column in a tabular input file')
    parser.add_argument('-o', '--output',
                        default=None,
                        help='Output file for query results')
    parser.add_argument('-O', '--output2',
                        default='iedb_results2',
                        help='Output file for secondary query results')
    parser.add_argument('-t', '--timeout',
                        type=int,
                        default=600,
                        help='Seconds to wait for server response')
    parser.add_argument('-r', '--retries',
                        type=int,
                        default=5,
                        help='Number of times to retry server query')
    parser.add_argument('-S', '--sleep',
                        type=int,
                        default=300,
                        help='Seconds to wait between retries')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False,
                        help='Turn on wrapper debugging to stderr')
    args = parser.parse_args()

    aapat = '^[ABCDEFGHIKLMNPQRSTVWY]+$'

    if args.prediction != 'bcell':
        if not args.allele and not args.allelefile:
            warn_err('-a allele or -A allelefile required\n', exit_code=1)

    if not (args.sequence or args.input):
        warn_err('NO Sequences given: '
                 + 'either -s sequence or -i input_file is required\n',
                 exit_code=1)

    if args.output is not None:
        try:
            outputPath = os.path.abspath(args.output)
            outputFile = open(outputPath, 'w')
        except Exception as e:
            warn_err("Unable to open output file: %s\n" % e, exit_code=1)
    else:
        outputFile = sys.stdout

    # params
    alleles = []
    lengths = []
    if args.prediction == 'bcell' and args.window_size is not None:
        lengths.append(str(args.window_size))
    else:
        if args.allelefile:
            (alleles, lengths) = parse_alleles(args.allelefile, args.length)
        if args.allele:
            for i, allele in enumerate(args.allele):
                alleles.append(allele)
                alen = args.length[i] if i < len(args.length)\
                    else args.length[-1]
                lengths.append(alen)
    allele = ','.join(alleles) if alleles else None
    length = ','.join(lengths) if lengths else None
    method = args.method
    proteasome = args.proteasome if args.prediction == 'processcing' else None
    url = 'https://tools-cluster-interface.iedb.org/tools_api/%s/' %\
        args.prediction

    # results
    results = {'prediction': {'header': None, 'entries': []},
               'detail': {'header': None, 'entries': []}}

    if args.sequence:
        for i, seq in enumerate(args.sequence):
            seqid = 'pep_%d' % i
            query(url, args.prediction, seq, allele, length, results,
                  seqid=seqid, method=method, proteasome=proteasome,
                  timeout=args.timeout, retries=args.retries,
                  sleep=args.sleep, debug=args.debug)
    if args.input:
        try:
            fh = open(args.input, 'r')
            if args.column:  # tabular
                col = int(args.column)
                idcol = int(args.id_column) if args.id_column else None
                for i, line in enumerate(fh):
                    fields = line.rstrip('\r\n').split('\t')
                    if len(fields) > col:
                        seq = re.sub('[_*]', '', fields[col].strip())
                        if re.match(aapat, seq):
                            if idcol is not None and idcol < len(fields):
                                seqid = fields[idcol]
                            else:
                                seqid = 'pep_%d' % i
                            query(url, args.prediction, seq, allele, length,
                                  results, seqid=seqid,
                                  method=method, proteasome=proteasome,
                                  timeout=args.timeout, retries=args.retries,
                                  sleep=args.sleep, debug=args.debug)
                        else:
                            warn_err('Line %d, Not a peptide: %s\n' % (i, seq),
                                     exit_code=None)
            else:  # fasta
                seqid = None
                seq = ''
                for i, line in enumerate(fh):
                    if line.startswith('>'):
                        if seqid and len(seq) > 0:
                            query(url, args.prediction, seq, allele, length,
                                  results, seqid=seqid,
                                  method=method, proteasome=proteasome,
                                  timeout=args.timeout, retries=args.retries,
                                  sleep=args.sleep, debug=args.debug)
                        seqid = line[1:].strip()
                        seq = ''
                    else:
                        seq += line.strip()
                if seqid and len(seq) > 0:
                    query(url, args.prediction, seq, allele, length,
                          results, seqid=seqid,
                          method=method, proteasome=proteasome,
                          timeout=args.timeout, retries=args.retries,
                          sleep=args.sleep, debug=args.debug)
            fh.close()
        except Exception as e:
            warn_err("Unable to open input file: %s\n" % e, exit_code=1)

    if results['prediction']['header']:
        outputFile.write(results['prediction']['header'])
    for line in results['prediction']['entries']:
        outputFile.write(line)
    if results['detail']['entries']:
        if args.output2:
            try:
                outPath = os.path.abspath(args.output2)
                outFile = open(outPath, 'w')
            except Exception as e:
                warn_err("Unable to open output file: %s\n" % e, exit_code=1)
        else:
            outFile = sys.stdout
        if results['detail']['header']:
            outFile.write(results['detail']['header'])
        for line in results['detail']['entries']:
            outFile.write(line)


if __name__ == "__main__":
    __main__()
