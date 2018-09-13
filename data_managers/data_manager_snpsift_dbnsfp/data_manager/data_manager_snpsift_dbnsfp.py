#!/usr/bin/env python

import gzip
import json
import optparse
import os
import os.path
import re
import shutil
import sys
import urllib
import zipfile

from pysam import ctabix

"""
# Install dbNSFP databases
# from DbNsfp site
  # Download dbNSFP database
    $ wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.4.zip
  # Uncompress
    $ unzip dbNSFP2.4.zip
  # Create a single file version
    $ (head -n 1 dbNSFP2.4_variant.chr1 ; cat dbNSFP2.4_variant.chr* | grep -v "^#") > dbNSFP2.4.txt
  # Compress using block-gzip algorithm
    bgzip dbNSFP2.4.txt
  # Create tabix index
    tabix -s 1 -b 2 -e 2 dbNSFP2.4.txt.gz

data_table:

    <table name="snpsift_dbnsfps" comment_char="#">
        <columns>key, build, name, value, annotations</columns>
        <file path="tool-data/snpsift_dbnsfps.loc" />
    </table>

#id     build   description     path    annotations
#GRCh37_dbNSFP2.4       GRCh37  GRCh37 dbNSFP2.4        /depot/snpeff/dbNSFP2.4.gz  SIFT_pred,Uniprot_acc
#GRCh38_dbNSFP2.7       GRCh38  GRCh38 dbNSFP2.7        /depot/snpeff/dbNSFP2.7.gz  SIFT_pred,Uniprot_acc
"""

data_table = 'snpsift_dbnsfps'
softgenetics_url = 'ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/'
dbNSFP_file_pat = '(dbNSFP(.*)_variant|dbscSNV(.*)).chr(.*)'
tokenize = re.compile(r'(\d+)|(\D+)').findall
dbNSFP_name_pat = 'dbNSFP(v|_light)?(\d*).*?'


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def get_nsfp_genome_version(name):
    genome_version = 'hg19'
    dbNSFP_name_pat = '(dbscSNV|dbNSFP(v|_light)?)(\d*).*?'
    m = re.match(dbNSFP_name_pat, name)
    if m:
        (base, mid, ver) = m.groups()
        if base == 'dbscSNV':
            genome_version = 'hg19'
        else:
            genome_version = 'hg38' if ver == '3' else 'hg19' if ver == '2' else 'hg18'
    return genome_version


def get_annotations(gzip_path):
    annotations = None
    fh = None
    try:
        fh = gzip.open(gzip_path, 'r')
        buf = fh.read(10000)
        lines = buf.splitlines()
        headers = lines[0].split('\t')
        annotations = ','.join([x.strip() for x in headers[4:]])
    except Exception as e:
        stop_err('Error Reading annotations %s : %s' % (gzip_path, e))
    finally:
        if fh:
            fh.close()
    return annotations


def tabix_file(input_fname, output_fname):
    print >> sys.stdout, "tabix_file: %s -> %s" % (input_fname, output_fname)
    ctabix.tabix_compress(input_fname, output_fname, force=True)
    # Column indices are 0-based.
    ctabix.tabix_index(output_fname, seq_col=0, start_col=1, end_col=1)


def natural_sortkey(string):
    return tuple(int(num) if num else alpha for num, alpha in tokenize(string))


def download_dbnsfp_database(url, output_file):
    dbnsfp_tsv = None
    file_path = 'downloaded_file'
    urllib.urlretrieve(url, file_path)
    with zipfile.ZipFile(file_path, 'r') as my_zip:
        dbnsfp_tsv = output_file if output_file else 'dbnsfp_tsv'
        wtr = open(dbnsfp_tsv, 'w')
        allfiles = [info.filename for info in my_zip.infolist()]
        files = [f for f in allfiles if re.match(dbNSFP_file_pat, f)]
        files = sorted(files, key=natural_sortkey)
        for j, file in enumerate(files):
            tempfiles = []
            tempfiles.append(file + "_%d" % len(tempfiles))
            tfh = open(tempfiles[-1], 'w')
            lastpos = None
            fh = my_zip.open(file, 'rU')
            for i, line in enumerate(fh):
                if i == 0:
                    if j == 0:
                        wtr.write(line)
                    continue
                else:
                    pos = int(line.split('\t')[1])
                    if lastpos and pos < lastpos:
                        tfh.close()
                        tempfiles.append(file + "_%d" % len(tempfiles))
                        tfh = open(tempfiles[-1], 'w')
                        print >> sys.stderr, "%s [%d] pos: %d < %d" % (file, i, pos, lastpos)
                    lastpos = pos
                tfh.write(line)
            tfh.close()
            if len(tempfiles) == 1:
                with open(tempfiles[0], 'r') as tfh:
                    wtr.writelines(tfh.readlines())
            else:
                tfha = [open(temp, 'r') for temp in tempfiles]
                lines = [tfh.readline() for tfh in tfha]
                curpos = [int(line.split('\t')[1]) for line in lines]
                while len(tfha) > 0:
                    k = curpos.index(min(curpos))
                    wtr.write(lines[k])
                    line = tfha[k].readline()
                    if line:
                        lines[k] = line
                        curpos[k] = int(line.split('\t')[1])
                    else:
                        tfha[k].close()
                        del tfha[k]
                        del lines[k]
                        del curpos[k]
    return dbnsfp_tsv


def main():
    # Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option('-g', '--dbkey', dest='dbkey', action='store', type="string", default=None, help='dbkey genome version')
    parser.add_option('-n', '--db_name', dest='db_name', action='store', type="string", default=None, help='A name for a history snpsiftdbnsfp dataset')
    parser.add_option('-s', '--softgenetics', dest='softgenetics', action='store', type="string", default=None, help='A name for softgenetics dbNSFP file')
    parser.add_option('-H', '--snpsiftdbnsfp', dest='snpsiftdbnsfp', action='store', type="string", default=None, help='A history snpsiftdbnsfp dataset')
    parser.add_option('-T', '--dbnsfp_tabular', dest='dbnsfp_tabular', action='store', type="string", default=None, help='A history dbnsfp_tabular dataset')
    (options, args) = parser.parse_args()

    filename = args[0]
    params = json.loads(open(filename).read())
    target_directory = params['output_data'][0]['extra_files_path']
    if not os.path.exists(target_directory):
        os.mkdir(target_directory)
    data_manager_dict = {}
    genome_version = options.dbkey if options.dbkey else 'unknown'
    dbnsfp_tsv = None
    db_name = None
    bzip_path = None
    if options.softgenetics:
        dbnsfp_url = softgenetics_url + options.softgenetics
        db_name = options.db_name if options.db_name else re.sub('\.zip$', '', options.softgenetics)
        genome_version = get_nsfp_genome_version(options.softgenetics)
        tsv = db_name + '.tsv'
        dbnsfp_tsv = download_dbnsfp_database(dbnsfp_url, tsv)
    elif options.dbnsfp_tabular:
        db_name = options.db_name
        dbnsfp_tsv = options.dbnsfp_tabular
    elif options.snpsiftdbnsfp:
        (dirpath, bgzip_name) = os.path.split(options.snpsiftdbnsfp)
        idxpath = options.snpsiftdbnsfp + '.tbi'
        shutil.copy(options.snpsiftdbnsfp, target_directory)
        shutil.copy(idxpath, target_directory)
        bzip_path = os.path.join(target_directory, bgzip_name)
        db_name = re.sub('(.txt)?.gz$', '', bgzip_name)
    else:
        stop_err('Either --softgenetics or --dbnsfp_tabular required')
    if dbnsfp_tsv:
        bgzip_name = '%s.txt.gz' % db_name
        bzip_path = os.path.join(target_directory, bgzip_name)
        tabix_file(dbnsfp_tsv, bzip_path)
    annotations = get_annotations(bzip_path)
    # Create the SnpSift dbNSFP Reference Data
    data_table_entry = dict(key='%s_%s' % (genome_version, db_name), build=genome_version, name='%s %s' % (genome_version, db_name), value=bgzip_name, annotations=annotations)
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table] = data_manager_dict['data_tables'].get(data_table, [])
    data_manager_dict['data_tables'][data_table].append(data_table_entry)

    # save info to json file
    open(filename, 'wb').write(json.dumps(data_manager_dict))


if __name__ == "__main__":
    main()
