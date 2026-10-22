"""
Helper script to check if AMAS input files are interleaved.
"""
import argparse
import re
import sys


def check_phylip_interleaved(filepath):
    """Check if PHYLIP file is interleaved."""
    with open(filepath, encoding='utf-8') as f:
        # First line is header: ntax nchar
        header = next(f).strip().split()
        ntax = int(header[0])

        for idx, line in enumerate(f, 1):
            if line.strip():
                if idx > ntax:
                    return True

        return False


def check_nexus_interleaved(filepath):
    """Check if NEXUS file is interleaved."""
    in_data_block = False
    in_matrix = False
    ntax = None
    seq_lines = 0

    with open(filepath, encoding='utf-8') as f:
        for line in f:
            content = line.strip().lower()

            if not content:
                continue

            if in_matrix:
                if content == 'end;':
                    return seq_lines != ntax if ntax else False

                if content != ';':
                    seq_lines += 1
                    if ntax and seq_lines > ntax:
                        return True
                continue

            if not in_data_block:
                if content.startswith('begin'):
                    words = content.split()
                    if len(words) > 1 and (
                            words[1].startswith('data')
                            or words[1].startswith('characters')):
                        in_data_block = True
                continue

            if content.startswith('dimensions') and ntax is None:
                match = re.search(r'ntax=(\d+)', content)
                if match:
                    ntax = int(match.group(1))

            elif content.startswith('format'):
                if re.search(r'\binterleave(?:;|=yes;?)?\b', content):
                    return True

            elif content.startswith('matrix'):
                in_matrix = True

    return False


def check_fasta_interleaved(filepath):
    """FASTA files are not interleaved."""
    return False


def main():
    parser = argparse.ArgumentParser(
        description='Check if AMAS input files are interleaved'
    )
    parser.add_argument('input_files', nargs='+', help='Input sequence files')
    parser.add_argument('--format', required=True,
                        choices=['fasta', 'phylip', 'nexus'],
                        help='Input format')

    args = parser.parse_args()

    interleaved_status = []
    for filepath in args.input_files:
        if args.format == 'phylip':
            is_interleaved = check_phylip_interleaved(filepath)
        elif args.format == 'nexus':
            is_interleaved = check_nexus_interleaved(filepath)
        else:
            is_interleaved = check_fasta_interleaved(filepath)

        interleaved_status.append(is_interleaved)

    interleaved_status = list(set(interleaved_status))
    if len(interleaved_status) > 1:
        raise Exception("Error: Input files are a mix of interleaved/sequential formats")

    if interleaved_status[0]:
        print(f"{args.format}-int")
    else:
        print(args.format)

    return 0


if __name__ == '__main__':
    sys.exit(main())
