"""Preserve MaxBin failures except its explicit insufficient-marker outcome."""

import gzip
import subprocess
import sys
import tarfile
from pathlib import Path


NO_BINS = (
    "Marker gene search reveals that the dataset cannot be binned "
    "(the medium of marker gene number <= 1). Program stop."
)


def main(args):
    result = subprocess.run(["run_MaxBin.pl", *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sys.stdout.buffer.write(result.stdout)
    sys.stdout.buffer.flush()
    if result.stderr:
        sys.stderr.buffer.write(result.stderr)
        sys.stderr.buffer.flush()
    if result.returncode != 255 or result.stdout.decode(errors="replace").splitlines()[-1:] != [NO_BINS]:
        return result.returncode
    if result.stderr and result.stderr.strip():
        return result.returncode

    prefix = args[args.index("-out") + 1]
    # A partial result is not the known, initial seed-search failure.
    if list(Path(prefix).parent.glob(f"{Path(prefix).name}.[0-9]*.fasta")):
        return result.returncode
    contig = args[args.index("-contig") + 1]
    min_length = int(args[args.index("-min_contig_length") + 1])
    with open(contig, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(contig, "rt") as source, open(f"{prefix}.noclass", "w") as unclassified, open(
        f"{prefix}.tooshort", "w"
    ) as too_short:
        record = []
        length = 0
        for line in source:
            if line.startswith(">") and record:
                (unclassified if length >= min_length else too_short).writelines(record)
                record, length = [], 0
            record.append(line)
            if not line.startswith(">"):
                length += len(line.strip())
        if record:
            (unclassified if length >= min_length else too_short).writelines(record)
    for suffix in ("summary", "marker", "marker.pdf"):
        Path(f"{prefix}.{suffix}").touch()
    with tarfile.open(f"{prefix}.marker_of_each_bin.tar.gz", "w:gz"):
        pass
    if "-reassembly" in args:
        directory = Path(f"{prefix}.reassem")
        directory.mkdir(exist_ok=True)
        for name in ("out.summary", "out.reads.noclass", "N50.txt"):
            (directory / name).touch()
    print("No bins recovered: insufficient marker genes; returning an empty bin collection.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
