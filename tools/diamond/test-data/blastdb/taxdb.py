#!/usr/bin/env python3
"""
build_taxdb_from_nodes.py

Create taxdb.btd and taxdb.bti (NCBI/BLAST/ISAM format) from a pruned nodes.dmp
and optional names.dmp found in the current directory.

Usage:
    python3 build_taxdb_from_nodes.py

Output:
    taxdb.btd
    taxdb.bti

Notes:
- Writes integers in BIG-ENDIAN (network order) as required by the ISAM/NCBI format.
- The btd records are written as:
    scientific_name<TAB>common_name<TAB>blast_name<TAB>superkingdom_code
  with no reliance on newlines for delimitation (offsets define length).
"""
import struct
import sys
from collections import defaultdict

NODES_FILE = "../ncbi_taxonomy/nodes.dmp"
NAMES_FILE = "../ncbi_taxonomy/names.dmp"   # optional
OUT_BTD = "taxdb.btd"
OUT_BTI = "taxdb.bti"

TAXDB_MAGIC = 0x8739


# -------------------------
# Helpers
# -------------------------
def read_nodes(nodes_path):
    """Return dicts: parent[taxid]=parent_taxid, rank[taxid]=rank"""
    parent = {}
    rank = {}
    with open(nodes_path, encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[0])
                parent_tax = int(parts[1])
            except ValueError:
                continue
            parent[taxid] = parent_tax
            rank[taxid] = parts[2]
    return parent, rank


def read_names(names_path):
    """Return dict: names[taxid] = {'scientific':..., 'common':..., 'blast':...}"""
    names = defaultdict(lambda: {"scientific": "", "common": "", "blast": ""})
    with open(names_path, encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            name_txt = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                names[taxid]["scientific"] = name_txt
            elif name_class == "common name":
                names[taxid]["common"] = name_txt
            elif name_class == "blast name":
                names[taxid]["blast"] = name_txt
    return names


def infer_superkingdom_code(taxid, parent, rank, sci_name_lookup):
    """
    Walk ancestors until rank == 'superkingdom', then map name to code:
    B (Bacteria), A (Archaea), E (Eukaryota), V (Viruses), U (Unknown)
    """
    seen = set()
    cur = taxid
    while True:
        if cur in seen:
            return "Unknown"
        seen.add(cur)
        r = rank.get(cur, "")
        if r == "domain":
            name = sci_name_lookup.get(cur, "").lower()
            if "bacteria" in name or "eubacteria" in name:
                return "Bacteria"
            if "archaea" in name:
                return "Archaea"
            if "eukaryota" in name or "eukaryota" in name or "eukary" in name:
                return "Eukaryota"
            if "virus" in name or "viruses" in name:
                return "Viruses"
            return "Unknown"
        if cur not in parent:
            return "Unknown"
        cur = parent[cur]


def infer_blast_name(taxid, parent, lookup):
    """
    """
    seen = set()
    cur = taxid
    while True:
        if cur in seen:
            return "Unknown"
        seen.add(cur)
        name = lookup.get(cur, "").lower()

        if name:
            return name
        if cur not in parent:
            return "Unknown"
        cur = parent[cur]


# -------------------------
# Main
# -------------------------
def main():
    # Read nodes.dmp
    try:
        parent, rank = read_nodes(NODES_FILE)
    except FileNotFoundError:
        print(f"Error: {NODES_FILE} not found in current directory.", file=sys.stderr)
        sys.exit(2)

    # Read names.dmp if present
    try:
        names = read_names(NAMES_FILE)
    except FileNotFoundError:
        names = defaultdict(lambda: {"scientific": "", "common": "", "blast": ""})
        print("Warning: names.dmp not found. scientific_name will be set to the taxid.", file=sys.stderr)

    # Determine the taxids to write:
    # use taxids present in nodes.dmp (pruned set)
    taxids = sorted(parent.keys())

    if len(taxids) == 0:
        print("No taxids found in nodes.dmp; nothing to do.", file=sys.stderr)
        sys.exit(0)

    # Build scientific-name lookup for superkingdom inference
    sci_lookup = {}
    for tid, rec in names.items():
        sci_lookup[tid] = rec.get("scientific", "")

    # Build blast-name lookup blast name inference
    bla_lookup = {}
    for tid, rec in names.items():
        bla_lookup[tid] = rec.get("blast", "")

    # Build btd records and offsets
    offsets = []
    btd_buf = bytearray()
    for tid in taxids:
        offsets.append(len(btd_buf))
        rec = names.get(tid, {"scientific": "", "common": "", "blast": ""})
        sci = rec.get("scientific", "")
        com = rec.get("common", "")

        if not sci:
            # fallback: use numeric taxid as scientific name (ensures non-empty)
            sci = str(tid)

        # infer superkingdom code from nodes.dmp and names if possible
        sk = infer_superkingdom_code(tid, parent, rank, sci_lookup)
        bla = infer_blast_name(tid, parent, bla_lookup)

        # exactly 4 fields, tab-separated; no trailing newline required
        record = f"{sci}\t{com}\t{bla}\t{sk}"
        btd_buf.extend(record.encode("utf-8"))

    end_offset = len(btd_buf)

    # Write taxdb.btd
    with open(OUT_BTD, "wb") as fh:
        fh.write(btd_buf)

    # Write taxdb.bti
    with open(OUT_BTI, "wb") as fh:
        # header: magic, count (number of real taxids), reserved[4]
        # IMPORTANT: write all integers BIG-ENDIAN (>I)
        fh.write(struct.pack(">I", TAXDB_MAGIC))
        fh.write(struct.pack(">I", len(taxids)))     # n (real entries only)
        fh.write(struct.pack(">IIII", 0, 0, 0, 0))   # reserved

        # index entries: (taxid, offset) pairs
        for tid, off in zip(taxids, offsets):
            fh.write(struct.pack(">I", int(tid)))
            fh.write(struct.pack(">I", int(off)))

        # # sentinel entry: taxid=0, offset=end_of_btd
        # fh.write(struct.pack(">I", 0))
        # fh.write(struct.pack(">I", end_offset))

    # Summary
    print(f"Wrote {OUT_BTD} ({end_offset} bytes)")
    print(f"Wrote {OUT_BTI} (header + {len(taxids)} entries)")
    print(f"Taxids written: {len(taxids)}")


if __name__ == "__main__":
    main()
