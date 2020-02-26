#!/usr/bin/env python

import argparse
import os.path
import re
import sys
import time

from urllib.error import HTTPError
from urllib.parse import urlencode, unquote
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
all_methods = set(mhci_methods + mhcii_methods +
                  mhcnp_methods + bcell_methods)
prediction_lengths = {'mhci': range(8, 16),
                      'mhcii': range(11, 31),
                      'processing': range(8, 15),
                      'mhcnp': range(8, 12),
                      'bcell': range(8, 16)}


def warn_err(msg, exit_code=1):
    sys.stderr.write(msg)
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
    parser.add_argument('-l', '--length',
                        action="append",
                        default=[],
                        help='lengths for which to make predictions, ' +
                             '1 per allele')
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=None,
                        help='window_size for bcell prediction')
    parser.add_argument('-i', '--input',
                        default=None,
                        help='Input file for peptide sequences ' +
                             '(fasta or tabular)')
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

    if not args.allele and args.prediction != 'bcell':
        warn_err('-a allele required\n', exit_code=1)

    if not (args.sequence or args.input):
        warn_err('NO Sequences given: ' +
                 'either -s sequence or -i input_file is required\n',
                 exit_code=1)

    if args.output is not None:
        try:
            outputPath = os.path.abspath(args.output)
            outputFile = open(outputPath, 'w')
        except Exception as e:
            warn_err("Unable to open output file: %s\n" % e, exit_code=1)
    else:
        outputFile = sys.stdout

    url = 'http://tools-cluster-interface.iedb.org/tools_api/%s/' %\
        args.prediction
    len_param = 'length' if args.prediction != 'bcell' else 'window_size'

    # TODO parse alleles from the args.alleles file
    alleles = ','.join(args.allele) if args.prediction != 'bcell' else None
    lengths = ','.join(args.length)
    if args.prediction == 'bcell':
        lengths = args.window_size
    method = args.method
    proteasome = args.proteasome if args.prediction == 'processcing' else None
    global header
    header = None
    results = []
    global header2
    header2 = None
    results2 = []

    sequence_text = []

    def add_seq(seqid, seq):
        sid = seqid if seqid else "peptide%d" % len(sequence_text)
        sequence_text.append(">%s\n%s" % (sid, seq))

    def query(url, seq, allele, length, seqid=None, method='recommended'):
        global header
        global header2
        params = dict()
        if method:
            params['method'] = method.encode()
        if proteasome:
            params['proteasome'] = proteasome.encode()
        params['sequence_text'] = seq.encode()
        if allele is not None:
            params['allele'] = allele.encode()
        if length is not None:
            params[len_param] = str(length).encode()
        req_data = urlencode(params)
        if args.debug:
            print('url %s %s' % (url, unquote(req_data)), file=sys.stderr)
        retries = max(0, args.retries) + 1
        for retry in range(1, retries):
            response = None
            try:
                response = urlopen(url, data=req_data.encode('utf-8'),
                                   timeout=args.timeout)
                if response and response.getcode() == 200:
                    data = [line.decode() for line in response.readlines()]
                    if args.debug:
                        print(data, file=sys.stderr)
                    rslts = results
                    for ln, line in enumerate(data):
                        if line.lower().find('invalid') >= 0:
                            msg = '%s %s\n%s' % (url, unquote(req_data),
                                                 ''.join(data))
                            warn_err(msg, exit_code=1)
                        if line.find('eptide') > 0:
                            header = "#%s%s" %\
                                    ("ID\t" if seqid else "", line)
                            if args.debug:
                                print(header, file=sys.stderr)
                            continue
                        elif method == 'Bepipred' and line.find('Residue') > 0:
                            header2 = "#%s%s" %\
                                    ("ID\t" if seqid else "", line)
                            if args.debug:
                                print(header2, file=sys.stderr)
                            rslts = results2
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
                code = None if retry < args.retries else e.code
                warn_err("%d of %d Error connecting to IEDB server %s\n" %
                         (retry, retries, e),
                         exit_code=code)
                time.sleep(args.sleep)
            except Exception as e:
                warn_err("Error connecting to IEDB server %s\n" % e,
                         exit_code=3)

    if args.sequence:
        for i, seq in enumerate(args.sequence):
            query(url, seq, alleles, lengths, seqid=None, method=method)
    if args.input:
        try:
            fh = open(args.input, 'r')
            if args.column:  # tabular
                col = int(args.column)
                idcol = int(args.id_column) if args.id_column else None
                for i, line in enumerate(fh):
                    fields = line.split('\t')
                    if len(fields) > col:
                        seq = re.sub('[_*]', '', fields[col])
                        if re.match(aapat, seq):
                            if idcol is not None and idcol < len(fields):
                                seqid = fields[idcol]
                            else:
                                seqid = None
                            query(url, seq, alleles, lengths,
                                  seqid=seqid, method=method)
                        else:
                            warn_err('Line %d, Not a peptide: %s\n' % (i, seq),
                                     exit_code=None)
            else:  # fasta
                seqid = None
                seq = ''
                for i, line in enumerate(fh):
                    if line.startswith('>'):
                        if seqid and len(seq) > 0:
                            query(url, seq, alleles, lengths,
                                  seqid=seqid, method=method)
                        seqid = line[1:].strip()
                        seq = ''
                    else:
                        seq += line.strip()
                if seqid and len(seq) > 0:
                    query(url, seq, alleles, lengths,
                          seqid=seqid, method=method)
            fh.close()
        except Exception as e:
            warn_err("Unable to open input file: %s\n" % e, exit_code=1)

    if header:
        outputFile.write(header)
    for line in results:
        outputFile.write(line)
    if results2:
        if args.output2:
            try:
                outPath = os.path.abspath(args.output2)
                outFile = open(outPath, 'w')
            except Exception as e:
                warn_err("Unable to open output file: %s\n" % e, exit_code=1)
        else:
            outFile = sys.stdout
        if header2:
            outFile.write(header2)
        for line in results2:
            outFile.write(line)


if __name__ == "__main__":
    __main__()
