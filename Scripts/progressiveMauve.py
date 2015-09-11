import subprocess
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run progressiveMauve')
    parser.add_argument('fasta_file', type=file, help='Fasta file')
    parser.add_argument('--output', help='Output file', default='output')
    args = parser.parse_args()

    cmd = ['progressiveMauve', '--output=%s.xmfa' % args.output,
           args.fasta_file.name]
    subprocess.check_call(cmd)
