import string
import sys

print("Validating --maas mapping file ...")
allowed = set(string.ascii_letters + string.digits + '._-')
for line in open(sys.argv[1], 'r'):
    if line.startswith('#'):
        continue
    for s in line.strip().split('\t'):
        if not set(s).issubset(allowed):
            print(f"Invalid line in mapping file: {line}")
            sys.exit(1)
