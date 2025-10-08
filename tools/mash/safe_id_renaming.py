import json
import os
import re
import sys


input_file = sys.argv[1]
input_id=sys.argv[2]

safe_file_name = re.sub(r'[^\w\-_\.]', '_', input_id)

safe_file_name = re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fa\.gz|\.fasta\.gz|\.fastq|\.fq|\.fa|\.fasta)$', '', safe_file_name, flags=re.IGNORECASE)

os.rename(input_file, safe_file_name)
