#!/usr/bin/env bash
set -euo pipefail

# Regenerate fixture JSON files for all 5 AlphaGenome Galaxy tools.
# Runs each tool against the real API with the same params as the Galaxy <test>
# sections, then converts the output to the fixture format.
#
# Usage:
#   export ALPHAGENOME_API_KEY=<your-key>
#   bash tools/alphagenome/generate_test_fixtures.sh

if [[ -z "${ALPHAGENOME_API_KEY:-}" ]]; then
    echo "ERROR: ALPHAGENOME_API_KEY is not set" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

PASS=0
FAIL=0
MAX_FIXTURE_ROWS=30

run_tool() {
    local name="$1"
    shift
    echo "--- Running $name ---"
    if python3 "$@" 2>"$TMPDIR/${name}.log"; then
        tail -3 "$TMPDIR/${name}.log"
        echo "    OK"
        PASS=$((PASS + 1))
    else
        tail -10 "$TMPDIR/${name}.log" >&2
        echo "    FAILED (see log above)" >&2
        FAIL=$((FAIL + 1))
        return 1
    fi
}

# ── 1. variant_effect (VCF output → fixture JSON) ──────────────────────

run_tool "variant_effect" \
    alphagenome_variant_effect.py \
    --input  test-data/test_input.vcf \
    --output "$TMPDIR/variant_effect.vcf" \
    --api-key "$ALPHAGENOME_API_KEY" \
    --output-types RNA_SEQ \
    --sequence-length 128KB \
    --max-variants 3 \
    --verbose

python3 -c '
import json, sys

# Reverse the INFO_FIELD_MAP from the tool script
FIELD_TO_TYPE = {
    "AG_RNA_LFC": "RNA_SEQ",
    "AG_ATAC_LFC": "ATAC",
    "AG_CAGE_LFC": "CAGE",
    "AG_DNASE_LFC": "DNASE",
    "AG_HISTONE_LFC": "CHIP_HISTONE",
    "AG_TF_LFC": "CHIP_TF",
    "AG_SPLICE_LFC": "SPLICE_SITES",
    "AG_SPLICEUSE_LFC": "SPLICE_SITE_USAGE",
    "AG_SPLICEJNC_LFC": "SPLICE_JUNCTIONS",
    "AG_CONTACT_LFC": "CONTACT_MAPS",
    "AG_PROCAP_LFC": "PROCAP",
}

variants = []
for line in open(sys.argv[1]):
    if line.startswith("#"):
        continue
    fields = line.strip().split("\t")
    chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
    info = fields[7]
    scores = {}
    for kv in info.split(";"):
        if "=" not in kv:
            continue
        key, val = kv.split("=", 1)
        if key in FIELD_TO_TYPE:
            scores[FIELD_TO_TYPE[key]] = round(float(val), 6)
    variants.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "scores": scores})

with open(sys.argv[2], "w") as f:
    json.dump({"variants": variants}, f, indent=2)
    f.write("\n")
print(f"  -> wrote {len(variants)} variants to {sys.argv[2]}")
' "$TMPDIR/variant_effect.vcf" \
  "test-data/fixture_variant_effect.json"

# ── Helper: TSV → fixture JSON ──────────────────────────────────────────

tsv_to_fixture() {
    local tsv_file="$1"
    local json_file="$2"
    python3 -c '
import csv, json, sys

max_rows = int(sys.argv[3])

with open(sys.argv[1]) as f:
    reader = csv.reader(f, delimiter="\t")
    columns = next(reader)
    rows = []
    total = 0
    for row in reader:
        total += 1
        if len(rows) >= max_rows:
            continue
        typed = []
        for val in row:
            try:
                typed.append(int(val))
            except ValueError:
                try:
                    typed.append(float(val))
                except ValueError:
                    typed.append(val)
        rows.append(typed)

with open(sys.argv[2], "w") as f:
    json.dump({"columns": columns, "rows": rows}, f, indent=2)
    f.write("\n")
if total > max_rows:
    print(f"  -> wrote {len(columns)} columns, {len(rows)}/{total} rows (capped) to {sys.argv[2]}")
else:
    print(f"  -> wrote {len(columns)} columns, {len(rows)} rows to {sys.argv[2]}")
' "$tsv_file" "$json_file" "$MAX_FIXTURE_ROWS"
}

# ── 2. variant_scorer (TSV output) ─────────────────────────────────────

run_tool "variant_scorer" \
    alphagenome_variant_scorer.py \
    --input  test-data/test_input.vcf \
    --output "$TMPDIR/variant_scorer.tsv" \
    --api-key "$ALPHAGENOME_API_KEY" \
    --scorers RNA_SEQ \
    --sequence-length 16KB \
    --max-variants 3 \
    --verbose

tsv_to_fixture "$TMPDIR/variant_scorer.tsv" \
    "test-data/fixture_variant_scorer.json"

# ── 3. ism_scanner (TSV output) ────────────────────────────────────────

run_tool "ism_scanner" \
    alphagenome_ism_scanner.py \
    --input  test-data/test_regions.bed \
    --output "$TMPDIR/ism_scanner.tsv" \
    --api-key "$ALPHAGENOME_API_KEY" \
    --scorers RNA_SEQ \
    --sequence-length 16KB \
    --max-regions 1 \
    --max-region-width 10 \
    --verbose

tsv_to_fixture "$TMPDIR/ism_scanner.tsv" \
    "test-data/fixture_ism_scanner.json"

# ── 4. interval_predictor (TSV output) ─────────────────────────────────

run_tool "interval_predictor" \
    alphagenome_interval_predictor.py \
    --input  test-data/test_intervals.bed \
    --output "$TMPDIR/interval_predictor.tsv" \
    --api-key "$ALPHAGENOME_API_KEY" \
    --output-types RNA_SEQ \
    --sequence-length 16KB \
    --max-intervals 3 \
    --output-mode summary \
    --verbose

tsv_to_fixture "$TMPDIR/interval_predictor.tsv" \
    "test-data/fixture_interval_predictor.json"

# ── 5. sequence_predictor (TSV output) ─────────────────────────────────

run_tool "sequence_predictor" \
    alphagenome_sequence_predictor.py \
    --input  test-data/test_sequences.fa \
    --output "$TMPDIR/sequence_predictor.tsv" \
    --api-key "$ALPHAGENOME_API_KEY" \
    --output-types RNA_SEQ \
    --sequence-length 16KB \
    --max-sequences 2 \
    --output-mode summary \
    --verbose

tsv_to_fixture "$TMPDIR/sequence_predictor.tsv" \
    "test-data/fixture_sequence_predictor.json"

# ── Summary ─────────────────────────────────────────────────────────────

echo ""
echo "=== Done: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
