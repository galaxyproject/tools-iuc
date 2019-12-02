#!/usr/bin/env python
"""
Runs RAxML on a sequence file.
For use with RAxML version 8.2.4
"""
import fnmatch
import glob
import optparse


def getint(name):
    basename = name.partition('RUN.')
    if basename[2] != '':
        num = basename[2]
        return int(num)


def __main__():
    # Parse the primary wrapper's command line options
    parser = optparse.OptionParser()
    # (-b)
    parser.add_option("--bootseed", action="store", type="int", dest="bootseed", help="Random number for non-parametric bootstrapping")
    # (-N/#)
    parser.add_option("--number_of_runs", action="store", type="int", dest="number_of_runs", default=1, help="Number of alternative runs")
    # (-q)
    parser.add_option("--multiple_model", action="store", type="string", dest="multiple_model", help="Multiple Model File")
    # (-x)
    parser.add_option("--rapid_bootstrap_random_seed", action="store", type="int", dest="rapid_bootstrap_random_seed", help="Rapid Boostrap Random Seed")

    (options, args) = parser.parse_args()

    # Multiple runs - concatenate
    if options.number_of_runs > 1:
        if options.bootseed is None and options.rapid_bootstrap_random_seed is None:
            runfiles = glob.glob('RAxML*RUN*')
            runfiles.sort(key=getint)
            # Logs
            with open('RAxML_log.galaxy', 'w') as outfile:
                for filename in runfiles:
                    if fnmatch.fnmatch(filename, 'RAxML_log.galaxy.RUN.*'):
                        with open(filename, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
            # Parsimony Trees
            with open('RAxML_parsimonyTree.galaxy', 'w') as outfile:
                for filename in runfiles:
                    if fnmatch.fnmatch(filename, 'RAxML_parsimonyTree.galaxy.RUN.*'):
                        with open(filename, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
            # Results
            with open('RAxML_result.galaxy', 'w') as outfile:
                for filename in runfiles:
                    if fnmatch.fnmatch(filename, 'RAxML_result.galaxy.RUN.*'):
                        with open(filename, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
    # Multiple Model Partition Files
    if options.multiple_model:
        files = glob.glob('RAxML_bestTree.galaxy.PARTITION.*')
        if len(files) > 0:
            files.sort(key=getint)
            # Best Tree Partitions
            with open('RAxML_bestTreePartitions.galaxy', 'w') as outfile:
                for filename in files:
                    if fnmatch.fnmatch(filename, 'RAxML_bestTree.galaxy.PARTITION.*'):
                        with open(filename, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
        else:
            with open('RAxML_bestTreePartitions.galaxy', 'w') as outfile:
                outfile.write("No partition files were produced.\n")

        # Result Partitions
        files = glob.glob('RAxML_result.galaxy.PARTITION.*')
        if len(files) > 0:
            files.sort(key=getint)
            with open('RAxML_resultPartitions.galaxy', 'w') as outfile:
                for filename in files:
                    if fnmatch.fnmatch(filename, 'RAxML_result.galaxy.PARTITION.*'):
                        with open(filename, 'r') as infile:
                            for line in infile:
                                outfile.write(line)
        else:
            with open('RAxML_resultPartitions.galaxy', 'w') as outfile:
                outfile.write("No partition files were produced.\n")

    # DEBUG options
    with open('RAxML_info.galaxy', 'a') as infof:
        infof.write('\nOM: CLI options DEBUG START:\n')
        infof.write(options.__repr__())
        infof.write('\nOM: CLI options DEBUG END\n')


if __name__ == "__main__":
    __main__()
