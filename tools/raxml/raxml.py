#!/usr/bin/env python
"""
Runs RAxML on a sequence file.
For use with RAxML version 8.2.4
"""
import fnmatch
import glob
import optparse
import os
import subprocess
import sys


def stop_err(msg):
    sys.stderr.write("%s\n" % msg)
    sys.exit()


def getint(name):
    basename = name.partition('RUN.')
    if basename[2] != '':
        num = basename[2]
        return int(num)


def __main__():
    usage = "usage: %prog -T <threads> -s <input> -n <output> -m <model> [optional arguments]"

    # Parse the primary wrapper's command line options
    parser = optparse.OptionParser(usage=usage)
    # raxml binary name, hardcoded in the xml file
    parser.add_option("--binary", action="store", type="string", dest="binary", help="Command to run")
    # (-a)
    parser.add_option("--weightfile", action="store", type="string", dest="weightfile", help="Column weight file")
    # (-A)
    parser.add_option("--secondary_structure_model", action="store", type="string", dest="secondary_structure_model", help="Secondary structure model")
    # (-b)
    parser.add_option("--bootseed", action="store", type="int", dest="bootseed", help="Bootstrap random number seed")
    # (-c)
    parser.add_option("--numofcats", action="store", type="int", dest="numofcats", help="Number of distinct rate categories")
    # (-d)
    parser.add_option("--search_complete_random_tree", action="store_true", dest="search_complete_random_tree", help="Search with a complete random starting tree")
    # (-D)
    parser.add_option("--ml_search_convergence", action="store_true", dest="ml_search_convergence", help="ML search onvergence criterion")
    # (-e)
    parser.add_option("--model_opt_precision", action="store", type="float", dest="model_opt_precision", help="Model Optimization Precision (-e)")
    # (-E)
    parser.add_option("--excludefile", action="store", type="string", dest="excludefile", help="Exclude File Name")
    # (-f)
    parser.add_option("--search_algorithm", action="store", type="string", dest="search_algorithm", help="Search Algorithm")
    # (-F)
    parser.add_option("--save_memory_cat_model", action="store_true", dest="save_memory_cat_model", help="Save memory under CAT and GTRGAMMA models")
    # (-g)
    parser.add_option("--groupingfile", action="store", type="string", dest="groupingfile", help="Grouping File Name")
    # (-G)
    parser.add_option("--enable_evol_heuristics", action="store_true", dest="enable_evol_heuristics", help="Enable evol algo heuristics")
    # (-i)
    parser.add_option("--initial_rearrangement_setting", action="store", type="int", dest="initial_rearrangement_setting", help="Initial Rearrangement Setting")
    # (-I)
    parser.add_option("--posterior_bootstopping_analysis", action="store", type="string", dest="posterior_bootstopping_analysis", help="Posterior bootstopping analysis")
    # (-J)
    parser.add_option("--majority_rule_consensus", action="store", type="string", dest="majority_rule_consensus", help="Majority rule consensus")
    # (-k)
    parser.add_option("--print_branch_lengths", action="store_true", dest="print_branch_lengths", help="Print branch lengths")
    # (-K)
    parser.add_option("--multistate_sub_model", action="store", type="string", dest="multistate_sub_model", help="Multistate substitution model")
    # (-m)
    parser.add_option("--model_type", action="store", type="string", dest="model_type", help="Model Type")
    parser.add_option("--base_model", action="store", type="string", dest="base_model", help="Base Model")
    parser.add_option("--aa_empirical_freq", action="store_true", dest="aa_empirical_freq", help="Use AA Empirical base frequences")
    parser.add_option("--aa_search_matrix", action="store", type="string", dest="aa_search_matrix", help="AA Search Matrix")
    # (-n)
    parser.add_option("--name", action="store", type="string", dest="name", help="Run Name")
    # (-N/#)
    parser.add_option("--number_of_runs", action="store", type="int", dest="number_of_runs", help="Number of alternative runs")
    parser.add_option("--number_of_runs_bootstop", action="store", type="string", dest="number_of_runs_bootstop", help="Number of alternative runs based on the bootstop criteria")
    # (-M)
    parser.add_option("--estimate_individual_branch_lengths", action="store_true", dest="estimate_individual_branch_lengths", help="Estimate individual branch lengths")
    # (-o)
    parser.add_option("--outgroup_name", action="store", type="string", dest="outgroup_name", help="Outgroup Name")
    # (-O)
    parser.add_option("--disable_undetermined_seq_check", action="store_true", dest="disable_undetermined_seq_check", help="Disable undetermined sequence check")
    # (-p)
    parser.add_option("--random_seed", action="store", type="int", dest="random_seed", help="Random Number Seed")
    # (-P)
    parser.add_option("--external_protein_model", action="store", type="string", dest="external_protein_model", help="External Protein Model")
    # (-q)
    parser.add_option("--multiple_model", action="store", type="string", dest="multiple_model", help="Multiple Model File")
    # (-r)
    parser.add_option("--constraint_file", action="store", type="string", dest="constraint_file", help="Constraint File")
    # (-R)
    parser.add_option("--bin_model_parameter_file", action="store", type="string", dest="bin_model_parameter_file", help="Constraint File")
    # (-s)
    parser.add_option("--source", action="store", type="string", dest="source", help="Input file")
    # (-S)
    parser.add_option("--secondary_structure_file", action="store", type="string", dest="secondary_structure_file", help="Secondary structure file")
    # (-t)
    parser.add_option("--starting_tree", action="store", type="string", dest="starting_tree", help="Starting Tree")
    # (-T)
    parser.add_option("--threads", action="store", type="int", dest="threads", help="Number of threads to use")
    # (-u)
    parser.add_option("--use_median_approximation", action="store_true", dest="use_median_approximation", help="Use median approximation")
    # (-U)
    parser.add_option("--save_memory_gappy_alignments", action="store_true", dest="save_memory_gappy_alignments", help="Save memory in large gapped alignments")
    # (-V)
    parser.add_option("--disable_rate_heterogeneity", action="store_true", dest="disable_rate_heterogeneity", help="Disable rate heterogeneity")
    # (-W)
    parser.add_option("--sliding_window_size", action="store", type="string", dest="sliding_window_size", help="Sliding window size")
    # (-x)
    parser.add_option("--rapid_bootstrap_random_seed", action="store", type="int", dest="rapid_bootstrap_random_seed", help="Rapid Boostrap Random Seed")
    # (-y)
    parser.add_option("--parsimony_starting_tree_only", action="store_true", dest="parsimony_starting_tree_only", help="Generate a parsimony starting tree only")
    # (-z)
    parser.add_option("--file_multiple_trees", action="store", type="string", dest="file_multiple_trees", help="Multiple Trees File")

    (options, args) = parser.parse_args()
    cmd = []

    # Required parameters
    binary = options.binary
    cmd.append(binary)
    # Threads
    if options.threads > 1:
        threads = "-T %d" % options.threads
        cmd.append(threads)
    # Source
    source = "-s %s" % options.source
    cmd.append(source)
    # Hardcode to "galaxy" first to simplify the output part of the wrapper
    # name = "-n %s" % options.name
    name = "-n galaxy"
    cmd.append(name)
    # Model
    model_type = options.model_type
    base_model = options.base_model
    aa_search_matrix = options.aa_search_matrix
    aa_empirical_freq = options.aa_empirical_freq
    if model_type == 'aminoacid':
        model = "-m %s%s" % (base_model, aa_search_matrix)
        if aa_empirical_freq:
            model = "-m %s%s%s" % (base_model, aa_search_matrix, 'F')
        # (-P)
        if options.external_protein_model:
            external_protein_model = "-P %s" % options.external_protein_model
            cmd.append(external_protein_model)
    else:
        model = "-m %s" % base_model
    cmd.append(model)
    if model == "GTRCAT":
        # (-c)
        if options.numofcats:
            numofcats = "-c %d" % options.numofcats
            cmd.append(numofcats)
    # Optional parameters
    if options.number_of_runs_bootstop:
        number_of_runs_bootstop = "-N %s" % options.number_of_runs_bootstop
        cmd.append(number_of_runs_bootstop)
    else:
        number_of_runs_bootstop = ''
    if options.number_of_runs:
        number_of_runs_opt = "-N %d" % options.number_of_runs
        cmd.append(number_of_runs_opt)
    else:
        number_of_runs_opt = 0
    # (-a)
    if options.weightfile:
        weightfile = "-a %s" % options.weightfile
        cmd.append(weightfile)
    # (-A)
    if options.secondary_structure_model:
        secondary_structure_model = "-A %s" % options.secondary_structure_model
        cmd.append(secondary_structure_model )
    # (-b)
    if options.bootseed:
        bootseed = "-b %d" % options.bootseed
        cmd.append(bootseed)
    else:
        bootseed = 0
    # -C - doesn't work in pthreads version, skipped
    if options.search_complete_random_tree:
        cmd.append("-d")
    if options.ml_search_convergence:
        cmd.append("-D" )
    if options.model_opt_precision:
        model_opt_precision = "-e %f" % options.model_opt_precision
        cmd.append(model_opt_precision)
    if options.excludefile:
        excludefile = "-E %s" % options.excludefile
        cmd.append(excludefile)
    if options.search_algorithm:
        search_algorithm = "-f %s" % options.search_algorithm
        cmd.append(search_algorithm)
    if options.save_memory_cat_model:
        cmd.append("-F")
    if options.groupingfile:
        groupingfile = "-g %s" % options.groupingfile
        cmd.append(groupingfile)
    if options.enable_evol_heuristics:
        enable_evol_heuristics = "-G %f" % options.enable_evol_heuristics
        cmd.append(enable_evol_heuristics )
    if options.initial_rearrangement_setting:
        initial_rearrangement_setting = "-i %s" % options.initial_rearrangement_setting
        cmd.append(initial_rearrangement_setting)
    if options.posterior_bootstopping_analysis:
        posterior_bootstopping_analysis = "-I %s" % options.posterior_bootstopping_analysis
        cmd.append(posterior_bootstopping_analysis)
    if options.majority_rule_consensus:
        majority_rule_consensus = "-J %s" % options.majority_rule_consensus
        cmd.append(majority_rule_consensus)
    if options.print_branch_lengths:
        cmd.append("-k")
    if options.multistate_sub_model:
        multistate_sub_model = "-K %s" % options.multistate_sub_model
        cmd.append(multistate_sub_model)
    if options.estimate_individual_branch_lengths:
        cmd.append("-M")
    if options.outgroup_name:
        outgroup_name = "-o %s" % options.outgroup_name
        cmd.append(outgroup_name)
    if options.disable_undetermined_seq_check:
        cmd.append("-O")
    if options.random_seed:
        random_seed = "-p %d" % options.random_seed
        cmd.append(random_seed)
    multiple_model = None
    if options.multiple_model:
        multiple_model = "-q %s" % options.multiple_model
        cmd.append(multiple_model)
    if options.constraint_file:
        constraint_file = "-r %s" % options.constraint_file
        cmd.append(constraint_file)
    if options.bin_model_parameter_file:
        bin_model_parameter_file_name = "RAxML_binaryModelParameters.galaxy"
        os.symlink(options.bin_model_parameter_file, bin_model_parameter_file_name )
        bin_model_parameter_file = "-R %s" % options.bin_model_parameter_file
        # Needs testing. Is the hardcoded name or the real path needed?
        cmd.append(bin_model_parameter_file)
    if options.secondary_structure_file:
        secondary_structure_file = "-S %s" % options.secondary_structure_file
        cmd.append(secondary_structure_file)
    if options.starting_tree:
        starting_tree = "-t %s" % options.starting_tree
        cmd.append(starting_tree)
    if options.use_median_approximation:
        cmd.append("-u")
    if options.save_memory_gappy_alignments:
        cmd.append("-U")
    if options.disable_rate_heterogeneity:
        cmd.append("-V")
    if options.sliding_window_size:
        sliding_window_size = "-W %d" % options.sliding_window_size
        cmd.append(sliding_window_size)
    if options.rapid_bootstrap_random_seed:
        rapid_bootstrap_random_seed = "-x %d" % options.rapid_bootstrap_random_seed
        cmd.append(rapid_bootstrap_random_seed)
    else:
        rapid_bootstrap_random_seed = 0
    if options.parsimony_starting_tree_only:
        cmd.append("-y")
    if options.file_multiple_trees:
        file_multiple_trees = "-z %s" % options.file_multiple_trees
        cmd.append(file_multiple_trees)

    print "cmd list: ", cmd, "\n"

    full_cmd = " ".join(cmd)
    print "Command string: %s" % full_cmd

    try:
        proc = subprocess.Popen(args=full_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as err:
        sys.stderr.write("Error invoking command: \n%s\n\n%s\n" % (cmd, err))
        sys.exit(1)
    stdout, stderr = proc.communicate()
    return_code = proc.returncode
    if return_code:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write("Return error code %i from command:\n" % return_code)
        sys.stderr.write("%s\n" % cmd)
    else:
        sys.stdout.write(stdout)
        sys.stdout.write(stderr)

    # Multiple runs - concatenate
    if number_of_runs_opt > 0:
        if (bootseed == 0) and (rapid_bootstrap_random_seed == 0 ):
            runfiles = glob.glob('RAxML*RUN*')
            runfiles.sort(key=getint)
        # Logs
            outfile = open('RAxML_log.galaxy', 'w')
            for filename in runfiles:
                if fnmatch.fnmatch(filename, 'RAxML_log.galaxy.RUN.*'):
                    infile = open(filename, 'r')
                    filename_line = "%s\n" % filename
                    outfile.write(filename_line)
                    for line in infile:
                        outfile.write(line)
                    infile.close()
            outfile.close()
        # Parsimony Trees
            outfile = open('RAxML_parsimonyTree.galaxy', 'w')
            for filename in runfiles:
                if fnmatch.fnmatch(filename, 'RAxML_parsimonyTree.galaxy.RUN.*'):
                    infile = open(filename, 'r')
                    filename_line = "%s\n" % filename
                    outfile.write(filename_line)
                    for line in infile:
                        outfile.write(line)
                    infile.close()
            outfile.close()
        # Results
            outfile = open('RAxML_result.galaxy', 'w')
            for filename in runfiles:
                if fnmatch.fnmatch(filename, 'RAxML_result.galaxy.RUN.*'):
                    infile = open(filename, 'r')
                    filename_line = "%s\n" % filename
                    outfile.write(filename_line)
                    for line in infile:
                        outfile.write(line)
                    infile.close()
            outfile.close()
    # Multiple Model Partition Files
    if multiple_model:
        files = glob.glob('RAxML_bestTree.galaxy.PARTITION.*')
        if len(files) > 0:
            files.sort(key=getint)
            outfile = open('RAxML_bestTreePartitions.galaxy', 'w')
            # Best Tree Partitions
            for filename in files:
                if fnmatch.fnmatch(filename, 'RAxML_bestTree.galaxy.PARTITION.*'):
                    infile = open(filename, 'r')
                    filename_line = "%s\n" % filename
                    outfile.write(filename_line)
                    for line in infile:
                        outfile.write(line)
                    infile.close()
            outfile.close()
        else:
            outfile = open('RAxML_bestTreePartitions.galaxy', 'w')
            outfile.write("No partition files were produced.\n")
            outfile.close()

        # Result Partitions
        files = glob.glob('RAxML_result.galaxy.PARTITION.*')
        if len(files) > 0:
            files.sort(key=getint)
            outfile = open('RAxML_resultPartitions.galaxy', 'w')
            for filename in files:
                if fnmatch.fnmatch(filename, 'RAxML_result.galaxy.PARTITION.*'):
                    infile = open(filename, 'r')
                    filename_line = "%s\n" % filename
                    outfile.write(filename_line)
                    for line in infile:
                        outfile.write(line)
                    infile.close()
            outfile.close()
        else:
            outfile = open('RAxML_resultPartitions.galaxy', 'w')
            outfile.write("No partition files were produced.\n")
            outfile.close()

    # DEBUG options
    infof = open('RAxML_info.galaxy', 'a')
    infof.write('\nOM: CLI options DEBUG START:\n')
    infof.write(options.__repr__())
    infof.write('\nOM: CLI options DEBUG END\n')

if __name__ == "__main__":
    __main__()
