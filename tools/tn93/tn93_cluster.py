import json
import argparse
import shutil
import shlex
import subprocess


def cluster_to_fasta (*args):
    with open (args[0], "r") as fh:
        cluster_json = json.load (fh)
        with open (args[1], "w") as fh2:
            for c in cluster_json:
                if len(args) == 3:
                    if args[2] in c['members']:
                        cc = c['centroid'].split ('\n')
                        cc[0] = ">" + args[2]
                        print ("\n".join (cc), file = fh2)
                        continue
                print (c['centroid'], file = fh2)

    return (os.path.getmtime(out_file), len (cluster_json))

def main(arguments):
    shutil.copy(arguments.input, os.path.join(os.getcwd(), 'reference_msa.fa'))
    while True: # '-f', '-o', ref_cluster, '-t', "%g" % ref_threshold, ref_msa_cpy
        command = 'tn93-cluster -o clusters.json -t %g reference_msa.fa' % (ref_cluster, ref_threshold, ref_msa_cpy)
        input_stamp = run_command (command, ref_cluster, "extract representative clusters at threshold %g" % ref_threshold)
        input_stamp, cluster_count = cluster_to_fasta (ref_cluster, _ref_seq_name)
        print ('Found %d clusters' % cluster_count)
        if cluster_count <= settings.max_clusters:
            shutil.copy (ref_msa_cpy, ref_msa)
            break
        else:
            ref_threshold += ref_step
    os.remove (ref_msa_cpy)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')
    parser.add_argument('--input', help='Input MSA', required=True, type=str)
    parser.add_argument('--threshold', help='Threshold', required=True, type=str )
    parser.add_argument('--ambigs', help='Handle ambigs', required=True, type=str)
    parser.add_argument('--cluster-type', help='Cluster type', required=True, type=str)
    parser.add_argument('--overlap', help='Overlap', required=True, type=str)
    parser.add_argument('--fraction', help='Fraction', required=True, type=str)
    parser.add_argument('--max-clusters', help='Max clusters', dest='max_clusters', required=True, type=int)
    parser.add_argument('--min-clusters', help='Min clusters', dest='min_clusters', required=True, type=int)
    arguments = parser.parse_args()
    exit(main(arguments))