"""
Helper script to check if AMAS input files are interleaved.
Returns simple boolean and format info.
"""
import sys
import argparse
import io
import re

def check_phylip_interleaved(filepath):
    """Check if PHYLIP file is interleaved."""
    with io.open(filepath, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    
    if len(lines) < 2:
        return False
    
    # First line is header: ntax nchar
    header = lines[0].split()
    if len(header) != 2:
        return False
    
    try:
        ntax = int(header[0])
    except:
        return False
    
    data_lines = lines[1:]
    
    # If we have more data lines than ntax, it's interleaved
    if len(data_lines) > ntax:
        return True
    
    return False

def check_nexus_interleaved(filepath):
    """Check if NEXUS file is interleaved."""
    with io.open(filepath, 'r') as f:
        content = f.read().lower()
    
    # Primary check: Look for explicit interleave keyword
    # Matches: 'interleave', 'interleave;', 'interleave=yes'
    # Not match: 'interleave=no'
    # NOTE: may be other ways interleave can be expressed in Format
    pattern = r'\binterleave(;|=yes)?\b'
    if re.search(pattern, content):
        return True
    
    # Backup check: Parse dimensions and count data lines
    # Look for dimensions: "Dimensions ntax=10 nchar=234"
    dim_match = re.search(r'dimensions\s+ntax=(\d+)\s+nchar=(\d+)', 
                          content)
    if not dim_match:
        return False
    
    ntax = int(dim_match.group(1))
    
    # Find Matrix section
    matrix_match = re.search(r'\bmatrix\b', content)
    if not matrix_match:
        return False
    
    matrix_start = matrix_match.end()
    
    # Find END; after matrix
    end_match = re.search(r'\bend\s*;', content[matrix_start:])
    if not end_match:
        return False
    
    matrix_end = matrix_start + end_match.start()
    
    matrix_content = content[matrix_start:matrix_end]
    
    lines = matrix_content.split('\n')
    data_line_count = 0
    for line in lines:
        stripped = line.strip()
        if stripped and stripped != ';':
            data_line_count += 1
    
    # If more data lines than ntax, it's interleaved
    return data_line_count > ntax

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
    
    # Check all files
    interleaved_status = []
    for filepath in args.input_files:
        if args.format == 'phylip':
            is_interleaved = check_phylip_interleaved(filepath)
        elif args.format == 'nexus':
            is_interleaved = check_nexus_interleaved(filepath)
        else:  # fasta
            is_interleaved = check_fasta_interleaved(filepath)
        
        interleaved_status.append(is_interleaved)
    
    # Determine overall status
    if all(interleaved_status):
        print("True")
        return 0
    else:
        print("False")
        return 0

if __name__ == '__main__':
    sys.exit(main())