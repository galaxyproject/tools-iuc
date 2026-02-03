import re
import sys


def main(in_file, bins_file, out_file, pattern):
    members = {}
    if bins_file != 'None':
        with open(bins_file) as bins:
            next(bins)
            for m in bins:
                name, binNr = m.split('\t')
                contig = name.strip()
                members[contig] = "bin_" + binNr.strip()

    with open(in_file, 'r') as f, open(out_file, 'w') as g:
        print(f"using pattern '{pattern}'")
        g.write("protein_id,contig_id,keywords\n")
        # Patterns: prodigal: /^>(.*?)_([0-9]*) #/       phanotate: /^>(.*?)_CDS_([0-9]*) /
        for line in f:
            if line.startswith(">"):
                match = re.match(pattern, line)
                if not match:
                    print("failed to match", line)
                protein = match.group(1)
                contig = match.group(2)
                if contig in members:
                    contig = members[contig]
                g.write(f"{protein},{contig},None\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
