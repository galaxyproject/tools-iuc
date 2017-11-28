#!/usr/bin/env python
"""
Input: DNA Fasta File
Output: Tabular
Return Tabular File with predicted ORF's
Bjoern Gruening
"""
import os
import shutil
import subprocess
import sys
import tempfile

from glimmer2seq import glimmer2seq


def main():
    genome_seq_file = sys.argv[1]
    outfile_classic_glimmer = sys.argv[2]
    outfile_ext_path = sys.argv[3]
    oufile_genes = sys.argv[8]

    tag = 'glimmer_non_knowlegde_based_prediction'
    tempdir = tempfile.gettempdir()

    trainingset = os.path.join(tempdir, tag + ".train")
    icm = os.path.join(tempdir, tag + ".icm")

    longorfs = tempfile.NamedTemporaryFile()
    trainingset = tempfile.NamedTemporaryFile()
    icm = tempfile.NamedTemporaryFile()

    # glimmeropts = "-o0 -g110 -t30 -l"
    glimmeropts = "-o%s -g%s -t%s" % (sys.argv[4], sys.argv[5], sys.argv[6])
    if sys.argv[7] == "true":
        glimmeropts += " -l"

    """
        1. Find long, non-overlapping orfs to use as a training set
    """
    subprocess.Popen(["long-orfs", "-n", "-t", "1.15",
                      genome_seq_file, "-"], stdout=longorfs,
                     stderr=subprocess.PIPE).communicate()

    """
        2. Extract the training sequences from the genome file
    """
    subprocess.Popen(["extract", "-t",
                      genome_seq_file, longorfs.name], stdout=trainingset,
                     stderr=subprocess.PIPE).communicate()

    """
        3. Build the icm from the training sequences
    """

    # the "-" parameter is used to redirect the output to stdout
    subprocess.Popen(["build-icm", "-r", "-"],
                     stdin=open(trainingset.name), stdout=icm,
                     stderr=subprocess.PIPE).communicate()

    """
        Run Glimmer3
    """
    subprocess.Popen(["glimmer3", glimmeropts,
                      genome_seq_file, icm.name, os.path.join(tempdir, tag)],
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    if outfile_classic_glimmer.strip() != 'None':
        shutil.copyfile(os.path.join(tempdir, tag + ".predict"), outfile_classic_glimmer)
    if outfile_ext_path.strip() != 'None':
        shutil.copyfile(os.path.join(tempdir, tag + ".detail"), outfile_ext_path)

    glimmer2seq(os.path.join(tempdir, tag + ".predict"), genome_seq_file, oufile_genes)


if __name__ == "__main__":
    main()
