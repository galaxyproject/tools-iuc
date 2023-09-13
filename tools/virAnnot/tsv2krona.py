#!/usr/bin/env python3


# Name: ecsv2krona
# Author: Marie Lefebvre - INRAE
# Aims: Create krona charts based on csv files from module blast2ecsv.py


import argparse
import csv
import logging as log
import os


def main():
    options = _set_options()
    _set_log_level(options.verbosity)
    files_list = _create_abundance(options)
    _create_cmd(options, files_list)


def _create_abundance(options):
    """
    extract values from csv files
    and create abundance files
    """
    log.info("Calculating abundance.")
    abundance = dict()
    for i in options.input:
        file_path = i[0]
        file_name = file_name = os.path.basename(i[0])
        current_id = os.path.splitext(file_name)[0]
        log.info(current_id)
        abundance[current_id] = {}
        with open(file_path, 'r') as current_file:
            log.debug("Reading " + file_path)
            csv_reader = csv.reader(current_file, delimiter='\t')
            line_count = 0
            for row in csv_reader:
                if line_count == 0:
                    # headers
                    line_count += 1
                else:
                    # no annotation
                    if len(row) == 16:
                        current_reads_nb = int(row[2])
                        if row[14] in abundance[current_id]:
                            # add reads
                            abundance[current_id][row[14]]["reads_nb"] = abundance[current_id][row[14]]["reads_nb"] + current_reads_nb
                            abundance[current_id][row[14]]["contigs_nb"] = abundance[current_id][row[14]]["contigs_nb"] + 1
                        else:
                            abundance[current_id][row[14]] = dict({"reads_nb": current_reads_nb, "contigs_nb": 1})
    if not os.path.exists(options.output):
        os.mkdir(options.output)
    if not os.path.exists(options.output + "/abundances"):
        os.mkdir(options.output + "/abundances")
    abundance_files = []
    # write abundance files
    for samp in abundance:
        reads_file = open(options.output + "/abundances/" + samp + "_reads.txt", "w+")
        abundance_files.append(options.output + "/abundances/" + samp + "_reads.txt")
        for taxo in abundance[samp]:
            reads_file.write(str(abundance[samp][taxo]["reads_nb"]))
            reads_file.write("\t")
            reads_file.write("\t".join(taxo.split(";")))
            reads_file.write("\n")
        reads_file.close()
        contigs_file = open(options.output + "/abundances/" + samp + "_contigs.txt", "w+")
        abundance_files.append(options.output + "/abundances/" + samp + "_contigs.txt")
        for taxo in abundance[samp]:
            contigs_file.write(str(abundance[samp][taxo]["contigs_nb"]))
            contigs_file.write("\t")
            contigs_file.write("\t".join(taxo.split(";")))
            contigs_file.write("\n")
        contigs_file.close()
        log.debug("Abundance file created " + options.output + "/abundances/" + samp + "_contigs.txt")
    return abundance_files


def _create_cmd(options, files_list):
    """
    Create one krona file for all abundances
    """
    cmd = "ktImportText "
    for file in files_list:
        cmd += file + " "
    cmd += "-o " + options.output + "/text.krona.html"
    log.info("Creating HTML.")
    log.debug(cmd)
    os.system(cmd)


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='TAB delimited file with results of blast from blast2ecsv',
                        action='append', required=True, dest='input', nargs='+')
    parser.add_argument('-o', '--out', help='The output directory', action='store', type=str, default='./krona', dest='output')
    parser.add_argument('-v', '--verbosity', help='Verbose level', action='store', type=int, choices=[1, 2, 3, 4], default=1)
    args = parser.parse_args()
    return args


def _set_log_level(verbosity):
    if verbosity == 1:
        log_format = '%(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.INFO, format=log_format)
    elif verbosity == 3:
        log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.DEBUG, format=log_format)


if __name__ == "__main__":
    main()
