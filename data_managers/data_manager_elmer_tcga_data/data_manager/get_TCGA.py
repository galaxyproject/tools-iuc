#!/usr/bin/env python
# Nitesh Turaga

import sys
import os
import shutil
import optparse
from rpy2.robjects.packages import importr
# from json import loads, dumps


# Import R packages required for
base = importr('base')
elmer = importr('ELMER')


# def cleanup_before_exit(tmp_dir):
    # """Remove tmp director"""
    # if tmp_dir and os.path.exsits(tmp_dir):
        # shutil.rmtree(tmp_dir)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def get_elmer_params(params):
    disease = params['param_dict']['cancer_type']
    Meth = params['param_dict']['methylation']
    RNA = params['param_dict']['rna']
    Clinic = params['param_dict']['clinical']
    return disease, Meth, RNA, Clinic #, RNAtype


def downlaod_tcga_data(data_manager_dict, params, disease, Meth, RNA, Clinic):
    """
    Fetch TCGA data
    """
    try:
        basedir = "Data"
        elmer.getTCGA(disease, Meth, RNA, Clinic, basedir)
        print " {} successfully downloaded in : {}".format(disease,basedir)
    except Exception, e:  # Catch FTP errors
        print str(e)
    # Delete raw data from TCGA, This is stored in"Data/Raw/"
    raw_data_dir = os.path.join(basedir,"Raw")
    if os.path.exsits(raw_data_dir):
        shutil.rmtree(raw_data_dir)
    # Add entry to data table
    if os.path.exists(basedir):
        files = os.listdir(os.path.join(basedir,disease))
        meth_file = disease+"_meth.rda"
        if meth_file in files:
            data_table_entry = dict(value = meth_file, name = disease, path=meth_file)
    data_table_entry = dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_name )
    _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    #  parse command line option
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()

    filename = args[0]
    params = loads(open(filename).read())



if __name__ == "__main__":
    main()
