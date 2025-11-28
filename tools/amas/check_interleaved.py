"""
Helper script to check if AMAS input files are interleaved.
"""
import argparse
import io
import re
import sys

def check_phylip_interleaved(filepath):
    """Check if PHYLIP file is interleaved."""
    seq_lines = 0

    with io.open(filepath, 'r') as f:
        for i, line in enumerate(f):
            # First line is header: ntax nchar
            if i==0:
                header = line[0].split()
                ntax = int(header[0])
            elif seq_lines > ntax:
                return True
            else:
                if line.strip():
                    seq_lines += 1
        
        return False

def check_nexus_interleaved(filepath):
    """Check if NEXUS file is interleaved."""
    block_flag = False
    seq_flag = False
    seq_lines = 0

    with io.open(filepath, 'r') as f:
        for line in f:
            content = line.lower().strip()

            if not content:
                continue

            if seq_flag:
                # Could shorten this to just `if content == 'end;':`
                if content == 'end;' and (seq_lines == ntax):
                    return False
                if content and content != ';':
                    seq_lines += 1
                if seq_lines > ntax:
                    return True

            split_content = content.split()

            if block_flag:
                if split_content[0] == 'dimensions':
                    # Will check for 'ntax=(int)'; position may not always be the same
                    ntax = int(re.search(r'ntax=(\d+)', content).group(1))

                elif split_content[0] == 'format':
                    # NOTE: may be other ways interleave can be expressed in Format
                    # Will match 'interleave', 'interleave;', 'interleave=yes', 'interleave=yes;'
                    if re.search(r'\binterleave(;|=yes|=yes;)?\b', content):
                        return True
                
                elif split_content[0] == 'matrix':
                    seq_flag = True
            
            # Need to make sure the block isn't a TAXA block that contains taxa info but no sequence info
            if split_content[0] == 'begin':
                if (split_content[1] == 'data;') or (split_content[1] == 'characters;'):
                    block_flag = True


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
    
    # NOTE: Do we need to check all files?
    if all(interleaved_status):
        return 0  # Exit code 0 = interleaved
    else:
        return 1  # Exit code 1 = sequential

if __name__ == '__main__':
    sys.exit(main())