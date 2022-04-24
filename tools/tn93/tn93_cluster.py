import argparse
import json
import os
import shlex
import subprocess
import sys


def cluster_to_fasta(json_file, fasta_file, reference_name=None):
    with open(json_file, "r") as fh:
        cluster_json = json.load(fh)
        with open(fasta_file, "w") as fh2:
            for c in cluster_json:
                if reference_name is not None:
                    if reference_name in c['members']:
                        cc = c['centroid'].split('\n')
                        cc[0] = ">" + reference_name
                        print("\n".join(cc), file=fh2)
                        continue
                print(c['centroid'], file=fh2)

    return(os.path.getmtime(fasta_file), len(cluster_json))


def run_command(command):
    proc = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    result = proc.returncode
    if result != 0:
        print('Command `%s` failed with exit code %s\n' % (command, result), file=sys.stderr)
        print('--------------------- STDOUT ---------------------')
        print(stdout.decode().replace('\\n', '\n'))
        print('------------------- END STDOUT -------------------')
        print('--------------------- STDERR ---------------------', file=sys.stderr)
        print(stderr.decode().replace('\\n', '\n'), file=sys.stderr)
        print('------------------- END STDERR -------------------', file=sys.stderr)
    return(int(result))


def main(arguments):
    threshold = arguments.threshold
    step = threshold * 0.25
    with open(arguments.reference) as fh:
        for line in fh:
            if line[0] == '>':
                _ref_seq_name = line[1:].split(' ')[0].strip()
                break
    while threshold <= 1:
        command = 'tn93-cluster -o clusters.json -t %g -a %s -c %s -m json -l %d -g %f %s' % (threshold, arguments.ambigs, arguments.cluster_type, arguments.overlap, arguments.fraction, arguments.input)
        return_code = run_command(command)
        if return_code != 0:
            return return_code
        input_stamp, cluster_count = cluster_to_fasta('clusters.json', 'clusters.fa', _ref_seq_name)
        if cluster_count <= arguments.cluster_count:
            break
        else:
            threshold += step
        print('Found %d clusters at threshold %f' % (cluster_count, threshold))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')
    parser.add_argument('--input', help='Input MSA', required=True, type=str)
    parser.add_argument('--reference', help='Reference sequence', required=True, type=str)
    parser.add_argument('--output', help='Input MSA', required=True, type=str)
    parser.add_argument('--threshold', help='Threshold', required=True, type=float)
    parser.add_argument('--ambigs', help='Handle ambigs', required=True, type=str)
    parser.add_argument('--cluster-type', help='Cluster type', required=True, type=str)
    parser.add_argument('--overlap', help='Overlap', required=True, type=int)
    parser.add_argument('--fraction', help='Fraction', required=True, type=float)
    parser.add_argument('--cluster-count', help='Max query', required=True, type=int)
    parser.add_argument('--compressed', help='File to write compressed clusters to', required=True, type=str)
    arguments = parser.parse_args()
    exit(main(arguments))
